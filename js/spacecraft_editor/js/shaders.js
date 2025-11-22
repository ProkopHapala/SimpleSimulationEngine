
const Shaders = {
    nodeVertex: `
        precision highp float;
        
        // Standard Three.js uniforms/attributes for RawShaderMaterial
        uniform mat4 modelViewMatrix;
        uniform mat4 projectionMatrix;
        attribute vec3 position;
        attribute vec3 normal;

        // Instancing attributes
        attribute float aNodeID; // Index into uPosTex (Float because attributes are floats by default unless specified)

        // Uniforms
        uniform sampler2D uPosTex;
        uniform vec2 uTexSize; // (width, height)
        uniform float uScale;

        varying vec3 vNormal;
        varying vec3 vColor;

        // Helper to fetch position from texture
        vec4 fetchPos(float id) {
            float u = (id + 0.5) / uTexSize.x;
            float v = 0.5 / uTexSize.y; // Assuming 1D texture (Nx1)
            return texture2D(uPosTex, vec2(u, v));
        }

        void main() {
            vec4 nodeData = fetchPos(aNodeID);
            vec3 center = nodeData.xyz;
            float type = nodeData.w;

            // Billboard logic (Sphere Impostor base)
            // For now, just a simple mesh instance
            vec3 transformed = position * uScale + center;

            vNormal = normal;
            
            // Color based on type (simple debug coloring)
            if (type == 0.0) vColor = vec3(1.0, 0.2, 0.2); // Red
            else if (type == 1.0) vColor = vec3(0.2, 1.0, 0.2); // Green
            else vColor = vec3(0.2, 0.2, 1.0); // Blue

            gl_Position = projectionMatrix * modelViewMatrix * vec4(transformed, 1.0);
        }
    `,
    nodeFragment: `
        precision highp float;
        
        varying vec3 vNormal;
        varying vec3 vColor;

        void main() {
            // Simple lighting
            vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
            float diff = max(dot(normalize(vNormal), lightDir), 0.0);
            vec3 color = vColor * (0.2 + 0.8 * diff);
            
            gl_FragColor = vec4(color, 1.0);
        }
    `,
    girderVertex: `
        precision highp float;

        // Standard Three.js uniforms/attributes
        uniform mat4 modelViewMatrix;
        uniform mat4 projectionMatrix;
        uniform mat4 viewMatrix; // Needed for billboard/alignment if used, or just standard
        attribute vec3 position;
        attribute vec3 normal;

        // Instancing attributes
        attribute vec2 aNodeIDs; // [idA, idB]
        
        // Uniforms
        uniform sampler2D uPosTex;
        uniform vec2 uTexSize;
        uniform float uThickness;

        varying vec3 vNormal;

        vec4 fetchPos(float id) {
            float u = (id + 0.5) / uTexSize.x;
            float v = 0.5 / uTexSize.y;
            return texture2D(uPosTex, vec2(u, v));
        }

        void main() {
            vec3 posA = fetchPos(aNodeIDs.x).xyz;
            vec3 posB = fetchPos(aNodeIDs.y).xyz;

            vec3 dir = posB - posA;
            float len = length(dir);
            vec3 axis = normalize(dir);

            vec3 center = (posA + posB) * 0.5;

            // Transform local cylinder vertex
            // Local Y is along the length
            // CylinderBufferGeometry is created with height 1, centered at 0
            vec3 localPos = position;
            localPos.x *= uThickness;
            localPos.z *= uThickness;
            localPos.y *= len; // Scale length

            // Rotate to align with axis
            vec3 up = vec3(0.0, 1.0, 0.0);
            vec3 zaxis = normalize(cross(up, axis));
            vec3 xaxis = cross(axis, zaxis);
            
            // Handle parallel case
            if (length(zaxis) < 0.001) {
                zaxis = vec3(1.0, 0.0, 0.0);
                xaxis = vec3(0.0, 0.0, 1.0);
            }
            
            mat3 rot = mat3(xaxis, axis, zaxis);
            vec3 worldPos = rot * localPos + center;

            vNormal = rot * normal;
            gl_Position = projectionMatrix * viewMatrix * vec4(worldPos, 1.0);
        }
    `,
    girderFragment: `
        precision highp float;
        
        varying vec3 vNormal;

        void main() {
            vec3 lightDir = normalize(vec3(0.5, 1.0, 0.5));
            float diff = max(dot(normalize(vNormal), lightDir), 0.0);
            vec3 color = vec3(0.7, 0.7, 0.8) * (0.3 + 0.7 * diff);
            gl_FragColor = vec4(color, 1.0);
        }
    `
};
