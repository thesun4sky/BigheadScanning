#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#define PI 3.14159265

class Vector3D{

public:
	float x,y,z;
public:
	Vector3D() {
	}
	Vector3D(float ex, float why, float zee) {
		x = ex, y = why, z = zee;
	}
	
	/*Vector3D(float ex = 0, float why = 0, float zee = 0) {
		x = ex, y = why, z = zee;
	}*/

	~Vector3D() {}

	float getMagnitude() {
		return sqrtf(x*x + y*y + z*z);
	}

	Vector3D operator*(float num) const {
		return Vector3D(x*num, y*num, z*num);
	}
	
	Vector3D operator+(float num) const {
		return Vector3D(x+num, y+num, z+num);
	}
	//∫§≈Õ µ°º¿
	Vector3D operator+(const Vector3D &vec) const {
		return Vector3D(x + vec.x, y + vec.y, z + vec.z);
	}
	//∫§≈Õ ª¨º¿
	Vector3D operator-(const Vector3D &vec) const {
		return Vector3D(x - vec.x, y - vec.y, z - vec.z);
	}
	//∫§≈Õ ¥Î¿‘
	Vector3D operator=(const Vector3D &vec)  {
		return Vector3D(vec.x, vec.y, vec.z);
	}
	//∫§≈Õ ¡§±‘»≠
	void normalizeVector3D(void) {
		float mag = sqrtf(x*x + y*y + z*z);
		x /= mag;
		y /= mag;
		z /= mag;
	}
	//∫§≈Õ¿« ≥ª¿˚
	float dotVector3D(const Vector3D &vec) const {
		return x*vec.x + y*vec.y + z*vec.z;
	}
	//∫§≈Õ¿« ø‹¿˚
	Vector3D crossVector3D(const Vector3D &vec) const {
		return Vector3D(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x);
	}
	Vector3D operator*(const Vector3D &vec) const {
		return Vector3D(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x);
	}
	//≥ª¿˚¿ª ¿ÃøÎ«ÿ µŒ ∫§≈Õ∞£¿« ∞¢
	float angleBetweenVector3Ds(Vector3D &vec) {
		return (acos(dotVector3D(vec) / (getMagnitude() * vec.getMagnitude())) * (180 / PI));
	}

};

#endif