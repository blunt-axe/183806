using UnityEngine;
using Unity.Collections;
using Unity.Mathematics;
using System.Collections;
using System.Threading.Tasks;

public class StructScript : MonoBehaviour {
	// 储存直线与物体表面相交的交点，向外的法向量，以及 cloth 的点
	public struct ConcatInfo {
		public float3 position;
		public float3 normal;
	}
	
	// 一个距离约束，储存距离，劲度系数（位于 [0, 1] 内），两个点的编号
	public struct DistanceConstraintInfo {
		public float dist, stiffness;
		public int u, v;
	};
	
	// 一个碰撞约束，储存 ConcatInfo 和点的编号
	public struct CollisionConstraintInfo {
		public ConcatInfo concat;
		public int u;
	};

	// 一个球体。对于一个绝对光滑的球体，不需要考虑它的角动量。
	public struct Sphere {
		public float3 position;
		public float3 velocity;
		public float radius;
		public float mass;

		public bool Contains(float3 point) {
			return radius - math.distance(point, position) >= 0;
		}

		// The direction **DOESN'T NEED** normalization. 
		// 注意，direction 必须是这一步指向上一步的向量，不能 normalize
		public ConcatInfo InternalRayIntersection(float3 point, float3 direction) {
			// if (this.Contains(point + direction)) {
				direction = math.normalize(point - position);
			// }

			direction = math.normalize(direction);
			point -= position;

			float b = math.dot(point, direction);
			float c = radius * radius - math.lengthsq(point);
			float t = -b + math.sqrt(b * b + c);

			float3 intersect = point + t * direction;
			float3 normal = math.normalize(intersect);
			return new ConcatInfo(){
				position = intersect + position,
				normal = normal
			};
		}
	}
}