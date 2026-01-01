attribute vec3 a_instance_pos;
attribute vec4 a_instance_quat;
varying vec3 vColor;

mat3 quat_to_mat3(vec4 q) {
    float x2 = q.x + q.x; float y2 = q.y + q.y; float z2 = q.z + q.z;
    float xx = q.x * x2;  float xy = q.x * y2;  float xz = q.x * z2;
    float yy = q.y * y2;  float yz = q.y * z2;  float zz = q.z * z2;
    float wx = q.w * x2;  float wy = q.w * y2;  float wz = q.w * z2;
    return mat3(
        1.0 - (yy + zz), xy + wz,          xz - wy,
        xy - wz,         1.0 - (xx + zz),  yz + wx,
        xz + wy,         yz - wx,          1.0 - (xx + yy)
    );
}

void main() {
    vColor = color;
    mat3 R = quat_to_mat3(a_instance_quat);
    vec3 transformed = R * position + a_instance_pos;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(transformed, 1.0);
}
