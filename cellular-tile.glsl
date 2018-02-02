// Source: https://www.shadertoy.com/view/MscSDB by Shane
// The cellular tile routine. Draw a few objects (four spheres, in this case) using a minumum
// blend at various 3D locations on a cubic tile. Make the tile wrappable by ensuring the 
// objects wrap around the edges. That's it.
//
// Believe it or not, you can get away with as few as three spheres. If you sum the total 
// instruction count here, you'll see that it's way, way lower than 2nd order 3D Voronoi.
// Not requiring a hash function provides the biggest benefit, but there is also less setup.
// 
// The result isn't perfect, but 3D cellular tiles can enable you to put a Voronoi looking 
// surface layer on a lot of 3D objects for little cost.
//
float drawSphere(in vec3 p){

    // Anything that wraps the domain will suffice, so any of the following will work.
    // p = cos(p*3.14159)*0.5;
    // p = abs(cos(p*3.14159)*0.5);
    // p = fract(p)-.5;
    // return dot(p, p);

    // Other metrics to try.

    p = abs(fract(p)-.5);
    return dot(p, vec3(.5));

    // p = abs(fract(p)-.5);
    // return max(max(p.x, p.y), p.z);

    // p = cos(p*3.14159)*0.5; 
    // p = abs(cos(p*3.14159)*0.5);
    // p = abs(fract(p)-.5);
    // return max(max(p.x - p.y, p.y - p.z), p.z - p.x);
    // return min(min(p.x - p.y, p.y - p.z), p.z - p.x);
}

// Faster (I'm assuming), more streamlined version. See the comments below for an expanded explanation.
// The function below is pretty quick also, and can be expanded to include more spheres. This one
// takes advantage of the fact that only four object need sorting. With three spheres, it'd be even
// better.
float cellTile(in vec3 p) {
    // Draw four overlapping objects (spheres, in this case) at various positions throughout the tile.
    vec4 v, d; 
    d.x = drawSphere(p - vec3(.81, .62, .53));
    p.xy = vec2(p.y-p.x, p.y + p.x)*.7071;
    d.y = drawSphere(p - vec3(.39, .2, .11));
    p.yz = vec2(p.z-p.y, p.z + p.y)*.7071;
    d.z = drawSphere(p - vec3(.62, .24, .06));
    p.xz = vec2(p.z-p.x, p.z + p.x)*.7071;
    d.w = drawSphere(p - vec3(.2, .82, .64));

    v.xy = min(d.xz, d.yw), v.z = min(max(d.x, d.y), max(d.z, d.w)), v.w = max(v.x, v.y); 

    d.x =  min(v.z, v.w) - min(v.x, v.y); // Maximum minus second order, for that beveled Voronoi look. Range [0, 1].
    //d.x =  min(v.x, v.y);
    return (d.x*2.66); // Normalize... roughly.
}

#pragma glslify: export(cellTile)
