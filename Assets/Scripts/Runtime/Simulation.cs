using UnityEngine;
using Unity.Collections;
using Unity.Mathematics;
using System.Collections;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;

using static StructScript;

public class ClothAndSphereController : MonoBehaviour
{
	// 必要的一些常数
	public const float sqrt2 = 1.41421356f;
	public float3 gAcceleration = new float3(0f, -9.8f, 0f);

	// 仿真的帧率
	public const int frameRate = 50;
	public const float dt = 0.02f;

	// 仿真求解过程中的参数
	public const float dampingConstant = 0.005f;
	public const int solverIterations = 50;
	public const float springStiffness = 0.5f;
	public float3 sprintStiffnessMultiplier = new float3(1f, 1f, 1f);
	public const float collisionStiffness = 1f;

	// 球体参数
	public const float sphereRadius = 1.5f;
	public const float sphereMass = 2.5f;

    // 布料参数
    public const int clothWidth = 40;
    public const int clothHeight = 40;
    public const float clothSpacing = 0.25f;
	public const float clothVertexMass = 1f;

	// 两个 GameObject
    private GameObject _sphere;
    private GameObject cloth;

	// 仿真中唯一球体的信息
	Sphere sphere;
	
	// 每次调用 Update() 时，clothMesh 为当前帧的布料，我们把它复制给 clothMeshTemp，
	// clothMeshTemp 会被用于渲染，与此同时直接开始对下一帧 clothMesh 的计算。
	// 注意：对 Mesh.vertices[i] 的赋值是无效的，只有 Mesh.vertices = Vector3[] 是有效的
    private Mesh clothMesh;
	private Mesh clothMeshTemp;

	// 用 NativeList 储存每个质点的信息。为啥不用 float3[]？因为 github 原项目用的 NativeList。
	private NativeList<float> masses;
	// Wass is the inversion of mass, LOL
	private NativeList<float> wasses;
	private NativeList<bool> isPinned;
	private NativeList<float3> positions;
	private NativeList<float3> velocities;
	private NativeList<float3> predictPositions;
	private NativeList<float3> predictPositionsTemp;

	// 布料内部的距离约束
	private NativeList<DistanceConstraintInfo> distanceConstraints;	

	// 碰撞约束
	private List<CollisionConstraintInfo> collisionConstraints;

	// NativeList 清空
	public void Dispose() {
		masses.Dispose();
		wasses.Dispose();
		isPinned.Dispose();
		positions.Dispose();
		velocities.Dispose();
		predictPositions.Dispose();
		predictPositionsTemp.Dispose();
		distanceConstraints.Dispose();
	}

	void OnDestroy() {
		Dispose();
		Destroy(_sphere);
		Destroy(cloth);
	}

	// 在 Awake 里更改仿真帧率
	void Awake() {
		// 想要改帧率，垂直同步一定要关掉（虽然不知道为啥。。）
		QualitySettings.vSyncCount = 0;
		Application.targetFrameRate = frameRate;
	}

    void Start()
    {
		// 写入 sphere 信息
		sphere = new Sphere() {
			position = new float3(3.5f, 1.5f, 2f),
			velocity = new float3(0f, 0f, 0f),
			radius = sphereRadius,
			mass = sphereMass
		};

        // Create the sphere
        _sphere = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        _sphere.transform.position = new Vector3(sphere.position[0], sphere.position[1], sphere.position[2]);
        _sphere.transform.localScale = Vector3.one * (2f * sphere.radius);

		// 设置材质的颜色、阴影
        Renderer sphereRenderer = _sphere.GetComponent<Renderer>();
        sphereRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.On; // Enable shadow casting
        sphereRenderer.receiveShadows = false; // Sphere usually doesn't receive shadows
		sphereRenderer.material.color = Color.red; // Make the sphere red

        // Create the cloth object
        cloth = new GameObject("Cloth");
        cloth.transform.position = new Vector3(0, 0, 0);

        // Create a MeshFilter and MeshRenderer for the cloth
        MeshFilter meshFilter = cloth.AddComponent<MeshFilter>();
        MeshRenderer meshRenderer = cloth.AddComponent<MeshRenderer>();
        clothMesh = new Mesh();
		clothMeshTemp = new Mesh();

		// **注意**这里用的是 clothMeshTemp
        meshFilter.mesh = clothMeshTemp;

        // Assign a basic material to render the cloth
        // meshRenderer.material = new Material(Shader.Find("Standard"));
        meshRenderer.material.color = new Color(0.18f, 0.54f, 0.34f);

        // Generate the grid mesh for the cloth
        GenerateClothMesh();

        // Enable shadows for the cloth
        meshRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off; // Cloth doesn't cast shadows
        meshRenderer.receiveShadows = true; // Cloth receives shadows

/*
        // Add a light source for shadows
        GameObject lightObj = new GameObject("Directional Light");
        Light light = lightObj.AddComponent<Light>();
        light.type = LightType.Directional;
        light.transform.rotation = Quaternion.Euler(50f, -30f, 0f); // Adjust the light angle
        light.shadows = LightShadows.Soft; // Enable soft shadows
		*/
    }

	// 每次调用 Update() 时，clothMesh 为当前帧的布料，我们把它复制给 clothMeshTemp，
	// clothMeshTemp 会被用于渲染，与此同时直接开始对下一帧 clothMesh 的计算。

	Task lastTask;

    async void Update()
    {
		if (lastTask != null) {
			await lastTask;
		}
        _sphere.transform.position = new Vector3(sphere.position[0], sphere.position[1], sphere.position[2]);
		clothMeshTemp.vertices = clothMesh.vertices;
		clothMeshTemp.normals = clothMesh.normals;

		// UpdateClothMesh() 是一个**异步**函数，可以认为它独立运行（大概吧）。
		// 理论上，只要它的运行时间不超过 20ms（即 1/50fps），the mesh would be updated。
        lastTask = UpdateClothMesh();
    }

	void InitializeClothLists() {
        int vertexCount = (clothWidth + 1) * (clothHeight + 1);

		masses = new NativeList<float>(vertexCount, Allocator.Persistent);
		wasses = new NativeList<float>(vertexCount, Allocator.Persistent);
		isPinned = new NativeList<bool>(vertexCount, Allocator.Persistent);
		positions = new NativeList<float3>(vertexCount, Allocator.Persistent);
		velocities = new NativeList<float3>(vertexCount, Allocator.Persistent);
		predictPositions = new NativeList<float3>(vertexCount, Allocator.Persistent);
		predictPositionsTemp = new NativeList<float3>(vertexCount, Allocator.Persistent);
		masses.Resize(vertexCount, NativeArrayOptions.ClearMemory);
		wasses.Resize(vertexCount, NativeArrayOptions.ClearMemory);
		isPinned.Resize(vertexCount, NativeArrayOptions.ClearMemory);
		positions.Resize(vertexCount, NativeArrayOptions.ClearMemory);
		velocities.Resize(vertexCount, NativeArrayOptions.ClearMemory);
		predictPositions.Resize(vertexCount, NativeArrayOptions.ClearMemory);
		predictPositionsTemp.Resize(vertexCount, NativeArrayOptions.ClearMemory);

		distanceConstraints = new NativeList<DistanceConstraintInfo>(6 * vertexCount, Allocator.Persistent);

		collisionConstraints = new List<CollisionConstraintInfo>(vertexCount);
	}

	void InitializeDistanceConstraints() {
		for (int y = 0; y < clothHeight; ++y) {
			for (int x = 0; x <= clothWidth; ++x) {
				distanceConstraints.Add(new DistanceConstraintInfo(){
					dist = clothSpacing,
					stiffness = springStiffness * sprintStiffnessMultiplier[0],
					u = x + y * (clothWidth + 1),
					v = x + (y + 1) * (clothWidth + 1)
				});
			}
		}

		for (int y = 0; y <= clothHeight; ++y) {
			for (int x = 0; x < clothWidth; ++x) {
				distanceConstraints.Add(new DistanceConstraintInfo(){
					dist = clothSpacing,
					stiffness = springStiffness * sprintStiffnessMultiplier[0],
					u = x + y * (clothWidth + 1),
					v = (x + 1) + y * (clothWidth + 1)
				});
			}
		}

		for (int y = 0; y <= clothHeight - 2; ++y) {
			for (int x = 0; x <= clothWidth; ++x) {
				distanceConstraints.Add(new DistanceConstraintInfo(){
					dist = 2f * clothSpacing,
					stiffness = springStiffness * sprintStiffnessMultiplier[1],
					u = x + y * (clothWidth + 1),
					v = x + (y + 2) * (clothWidth + 1)
				});
			}
		}

		for (int y = 0; y <= clothHeight; ++y) {
			for (int x = 0; x <= clothWidth - 2; ++x) {
				distanceConstraints.Add(new DistanceConstraintInfo(){
					dist = 2f * clothSpacing,
					stiffness = springStiffness * sprintStiffnessMultiplier[1],
					u = x + y * (clothWidth + 1),
					v = (x + 2) + y * (clothWidth + 1)
				});
			}
		}

		for (int y = 0; y < clothHeight; ++y) {
			for (int x = 0; x < clothWidth; ++x) {
				distanceConstraints.Add(new DistanceConstraintInfo(){
					dist = sqrt2 * clothSpacing,
					stiffness = springStiffness * sprintStiffnessMultiplier[2],
					u = x + y * (clothWidth + 1),
					v = (x + 1) + (y + 1) * (clothWidth + 1)
				});
				distanceConstraints.Add(new DistanceConstraintInfo(){
					dist = sqrt2 * clothSpacing,
					stiffness = springStiffness,
					u = x + (y + 1) * (clothWidth + 1),
					v = (x + 1) + y * (clothWidth + 1)
				});
			}
		}
	}

    void GenerateClothMesh()
    {
        int vertexCount = (clothWidth + 1) * (clothHeight + 1);
        Vector3[] vertices = new Vector3[vertexCount];
        int[] triangles = new int[clothWidth * clothHeight * 6];
        Vector2[] uvs = new Vector2[vertexCount];

		InitializeClothLists();

        // Create vertices and UVs
        for (int y = 0; y <= clothHeight; y++)
        {
            for (int x = 0; x <= clothWidth; x++)
            {
                int index = x + y * (clothWidth + 1);
                vertices[index] = new Vector3(x * clothSpacing, 0, y * clothSpacing);
                uvs[index] = new Vector2((float)x / clothWidth, (float)y / clothHeight);

				masses[index] = clothVertexMass;
				wasses[index] = 1.0f / clothVertexMass;
				positions[index] = new float3(x * clothSpacing, 0, y * clothSpacing);
				velocities[index] = new float3(0f, 0f, 0f);

				// isPinned[index] = (x == 0 || x == clothWidth || y == 0 || y == clothHeight) && (x % 5 == 0 || x == clothWidth) && (y % 5 == 0 || y == clothHeight);
				isPinned[index] = (x == 0 || x == clothWidth) && (y == 0 || y == clothHeight);
            }
        }

		InitializeDistanceConstraints();

        // Create triangles
        int triIndex = 0;
        for (int y = 0; y < clothHeight; y++)
        {
            for (int x = 0; x < clothWidth; x++)
            {
                int topLeft = x + y * (clothWidth + 1);
                int topRight = topLeft + 1;
                int bottomLeft = x + (y + 1) * (clothWidth + 1);
                int bottomRight = bottomLeft + 1;

                // First triangle
                triangles[triIndex++] = topLeft;
                triangles[triIndex++] = bottomLeft;
                triangles[triIndex++] = topRight;

                // Second triangle
                triangles[triIndex++] = topRight;
                triangles[triIndex++] = bottomLeft;
                triangles[triIndex++] = bottomRight;
            }
        }

        // Assign mesh data
        clothMesh.vertices = vertices;
        clothMesh.triangles = triangles;
        clothMesh.uv = uvs;
        clothMesh.RecalculateNormals();

		// Assign mesh data for clothMeshTemp（写得有点丑，凑合看）
        Vector3[] _vertices = new Vector3[vertexCount];
        int[] _triangles = new int[clothWidth * clothHeight * 6];
        Vector2[] _uvs = new Vector2[vertexCount];
		for (int i = 0; i < vertexCount; ++i) {
			_vertices[i] = vertices[i];
			_uvs[i] = uvs[i];
		}
		for (int i = 0; i < triIndex; ++i) {
			_triangles[i] = triangles[i];
		}
        clothMeshTemp.vertices = _vertices;
        clothMeshTemp.triangles = _triangles;
        clothMeshTemp.uv = _uvs;
        clothMeshTemp.RecalculateNormals();
    }

    public async Task UpdateClothMesh()
    {
		sphere.velocity = (1 - dampingConstant) * (sphere.velocity + dt * gAcceleration);
		Sphere predictSphere = new Sphere() {
			position = sphere.position + dt * sphere.velocity,
			velocity = sphere.velocity,
			radius = sphere.radius,
			mass = sphere.mass
		};

		int vertexCount = masses.Length;
		int distanceConstraintsCount = distanceConstraints.Length;

		for (int i = 0; i < vertexCount; ++i) {
			velocities[i] = (1 - dampingConstant) * (velocities[i] + dt * gAcceleration);
			if (isPinned[i]) {
				velocities[i] = new float3(0f, 0f, 0f);
			}
			predictPositions[i] = positions[i] + dt * velocities[i];
		}

		// 添加碰撞约束
		collisionConstraints.Clear();
		for (int i = 0; i < vertexCount; ++i) {
			if (!isPinned[i] && predictSphere.Contains(predictPositions[i])) {
				collisionConstraints.Add(new CollisionConstraintInfo(){
					concat = predictSphere.InternalRayIntersection(predictPositions[i], dt * (predictSphere.velocity - velocities[i])),
					u = i
				});
			}
		}

		int collisionConstraintsCount = collisionConstraints.Count;	

		for (int i = 0; i < collisionConstraintsCount; ++i) {
			ConcatInfo concat = collisionConstraints[i].concat;
			// Debug.Log("Contains " + collisionConstraints[i].u + "Position: " + concat.position[0] + " " + concat.position[1] + " " + concat.position[2] + " Normal: " + concat.normal[0] + " " + concat.normal[1] + " " + concat.normal[2]);
		}
		
		// 进行迭代
		for (int _ = 0; _ < solverIterations; ++_) {
			for (int i = 0; i < vertexCount; ++i) {
				predictPositionsTemp[i] = predictPositions[i];
			}
			for (int i = 0; i < distanceConstraintsCount; ++i) {
				int u = distanceConstraints[i].u;
				int v = distanceConstraints[i].v;
				float3 p1 = predictPositions[u];
				float3 p2 = predictPositions[v];
				// Magically, p12 = p1 - p2 (LOL)
				float3 p12 = p1 - p2;
				float realDist = math.length(p12);
				float expectedDist = distanceConstraints[i].dist;
				float stiffness = distanceConstraints[i].stiffness;
				float3 deltaP1 = -(wasses[u] / (wasses[u] + wasses[v])) * (realDist - expectedDist) / realDist * p12;
				float3 deltaP2 = +(wasses[v] / (wasses[u] + wasses[v])) * (realDist - expectedDist) / realDist * p12;
				predictPositionsTemp[u] += stiffness * deltaP1;
				predictPositionsTemp[v] += stiffness * deltaP2;
			}
			for (int i = 0; i < collisionConstraintsCount; ++i) {
				int u = collisionConstraints[i].u;
				float3 position = collisionConstraints[i].concat.position;
				float3 normal = collisionConstraints[i].concat.normal;
				if (math.dot(normal, predictPositions[u] - position) < 0) {
					predictPositionsTemp[u] += -math.dot(normal, predictPositions[u] - position) * collisionStiffness * normal;
				}
			}
			for (int i = 0; i < vertexCount; ++i) {
				if (!isPinned[i]) {
					predictPositions[i] += (predictPositionsTemp[i] - predictPositions[i]) / 10f;
				}
			}
		}

		for (int i = 0; i < collisionConstraintsCount; ++i) {
			int u = collisionConstraints[i].u;
			float3 normal = collisionConstraints[i].concat.normal;
			float3 vec = positions[u] + dt * velocities[u] - predictPositions[u];
			vec = math.dot(normal, vec) * normal;
			predictSphere.velocity += masses[u] / predictSphere.mass * vec;
		}
		
		predictSphere.position = sphere.position + dt * (sphere.velocity + predictSphere.velocity) / 2;

		sphere = predictSphere;

		Vector3[] vertices = new Vector3[vertexCount];
		for (int i = 0; i < vertexCount; ++i) {
			velocities[i] = (predictPositions[i] - positions[i]) / dt;
			positions[i] = predictPositions[i];
			for (int j = 0; j < 3; ++j) {
				vertices[i][j] = positions[i][j];
			}
		}
		
		// Update clothMesh
		clothMesh.vertices = vertices;
        clothMesh.RecalculateNormals();

		// Thread.Sleep(500);

		await Task.Yield();
    }
}
