// Equivalent to acos(1/sqrt(3))
#define FACE_VERTEX_EDGE 0.955317

#pragma glslify: rotationMatrix = require(../rotation-matrix3)

// Plane technique
// source: https://www.shadertoy.com/view/MtV3Dy
float  plane(vec3 p, vec3 origin, vec3 normal){ 
   return dot(p - origin,normal);
}

float planeTetrahedron (vec3 p, float d) {
  const float dn =1.0/sqrt(3.0);
  p *= rotationMatrix(normalize(vec3(-1.0, 0.0, 1.0)), FACE_VERTEX_EDGE);

  //The tetrahedran is the intersection of four planes:
  float sd1 = plane(p,vec3(d,d,d) ,vec3(-dn,dn,dn)) ; 
  float sd2 = plane(p,vec3(d,-d,-d) ,vec3(dn,-dn,dn)) ;
  float sd3 = plane(p,vec3(-d,d,-d) ,vec3(dn,dn,-dn)) ;
  float sd4 = plane(p,vec3(-d,-d,d) ,vec3(-dn,-dn,-dn)) ;

  //max intersects shapes
  return max(max(sd1,sd2),max(sd3,sd4));
}

#pragma glslify: export(planeTetrahedron)
