import * as THREE from 'three';

export class Physics {
    constructor(gl, nParticles, bodyDef) {
        this.gl = gl;
        this.nParticles = nParticles;
        this.bodyDef = bodyDef;
        this.size = Math.ceil(Math.sqrt(nParticles));
        this.pingPong = 0;
        
        this.initShaders();
        this.initTextures();
        this.initFramebuffers();
        this.initGeometry();
    }

    async initShaders() {
        const response = await fetch('RigidBody.glsl');
        const fsSource = await response.text();
        const vsSource = `#version 300 es
            in vec2 position;
            void main() { gl_Position = vec4(position, 0, 1); }`;

        this.program = this.createProgram(vsSource, fsSource);
        this.locations = {
            u_tex_pos: this.gl.getUniformLocation(this.program, 'u_tex_pos'),
            u_tex_vel: this.gl.getUniformLocation(this.program, 'u_tex_vel'),
            u_tex_quat: this.gl.getUniformLocation(this.program, 'u_tex_quat'),
            u_tex_ang_vel: this.gl.getUniformLocation(this.program, 'u_tex_ang_vel'),
            u_tex_force: this.gl.getUniformLocation(this.program, 'u_tex_force'),
            u_tex_torque: this.gl.getUniformLocation(this.program, 'u_tex_torque'),
            u_dt: this.gl.getUniformLocation(this.program, 'u_dt'),
            u_gravity: this.gl.getUniformLocation(this.program, 'u_gravity'),
            u_point_count: this.gl.getUniformLocation(this.program, 'u_point_count'),
            u_points: this.gl.getUniformLocation(this.program, 'u_points'),
            u_inertia_inv: this.gl.getUniformLocation(this.program, 'u_inertia_inv')
        };
    }

    createProgram(vsSource, fsSource) {
        const gl = this.gl;
        const vs = gl.createShader(gl.VERTEX_SHADER);
        gl.shaderSource(vs, vsSource);
        gl.compileShader(vs);
        if (!gl.getShaderParameter(vs, gl.COMPILE_STATUS)) {
            console.error('VS Error:', gl.getShaderInfoLog(vs));
        }

        const fs = gl.createShader(gl.FRAGMENT_SHADER);
        gl.shaderSource(fs, fsSource);
        gl.compileShader(fs);
        if (!gl.getShaderParameter(fs, gl.COMPILE_STATUS)) {
            console.error('FS Error:', gl.getShaderInfoLog(fs));
        }

        const program = gl.createProgram();
        gl.attachShader(program, vs);
        gl.attachShader(program, fs);
        gl.linkProgram(program);
        if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
            console.error('Link Error:', gl.getProgramInfoLog(program));
        }
        return program;
    }

    initTextures() {
        this.textures = [
            this.createStateTextures(),
            this.createStateTextures()
        ];
        this.reset(0);
    }

    createStateTextures() {
        return {
            pos: this.createFloatTexture(),
            vel: this.createFloatTexture(),
            quat: this.createFloatTexture(),
            angVel: this.createFloatTexture(),
            force: this.createFloatTexture(),
            torque: this.createFloatTexture()
        };
    }

    createFloatTexture() {
        const gl = this.gl;
        const tex = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, tex);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, this.size, this.size, 0, gl.RGBA, gl.FLOAT, null);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        return tex;
    }

    initFramebuffers() {
        const gl = this.gl;
        this.fbos = [gl.createFramebuffer(), gl.createFramebuffer()];
        this.setupFBO(this.fbos[0], this.textures[0]);
        this.setupFBO(this.fbos[1], this.textures[1]);
    }

    setupFBO(fbo, texs) {
        const gl = this.gl;
        gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texs.pos, 0);
        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, texs.vel, 0);
        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT2, gl.TEXTURE_2D, texs.quat, 0);
        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT3, gl.TEXTURE_2D, texs.angVel, 0);
        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT4, gl.TEXTURE_2D, texs.force, 0);
        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT5, gl.TEXTURE_2D, texs.torque, 0);
        gl.drawBuffers([gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1, gl.COLOR_ATTACHMENT2, gl.COLOR_ATTACHMENT3, gl.COLOR_ATTACHMENT4, gl.COLOR_ATTACHMENT5]);
    }

    initGeometry() {
        const gl = this.gl;
        this.quadVAO = gl.createVertexArray();
        gl.bindVertexArray(this.quadVAO);
        const posBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, posBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1,-1, 1,-1, -1,1, -1,1, 1,-1, 1,1]), gl.STATIC_DRAW);
        gl.enableVertexAttribArray(0);
        gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    }

    reset(scatterRad = 0, speedBase = 1.0, speedSpread = 1.0) {
        const gl = this.gl;
        const count = this.size * this.size;
        const pos  = new Float32Array(count * 4);
        const quat = new Float32Array(count * 4);
        const vel  = new Float32Array(count * 4);
        const zero = new Float32Array(count * 4);

        const rotateVecByQuat = (v, q) => {
            const x = q[0], y = q[1], z = q[2], w = q[3];
            const uvx  = 2 * (y * v[2] - z * v[1]);
            const uvy  = 2 * (z * v[0] - x * v[2]);
            const uvz  = 2 * (x * v[1] - y * v[0]);
            const uuvx = 2 * (y * uvz - z * uvy);
            const uuvy = 2 * (z * uvx - x * uvz);
            const uuvz = 2 * (x * uvy - y * uvx);
            return [v[0] + w * uvx + uuvx, v[1] + w * uvy + uuvy, v[2] + w * uvz + uuvz];
        };

        for (let i = 0; i < this.nParticles; i++) {
            // FIXED SEED FOR PARTICLE 0
            if (i === 0) {
                pos[0] = 0; pos[1] = 5.0; pos[2] = 0; pos[3] = 1.0;
                quat[0] = 0; quat[1] = 0; quat[2] = 0; quat[3] = 1.0;
                vel[0] = 0; vel[1] = 0; vel[2] = speedBase; vel[3] = 0;
            } else {
                pos[i*4  ] =       (Math.random() - 0.5) * 5.0;
                pos[i*4+1] = 5.0 + (Math.random() - 0.5) * 5.0;
                pos[i*4+2] =       (Math.random() - 0.5) * 5.0;
                pos[i*4+3] = 1.0;

                let axis = [Math.random()-0.5, Math.random()-0.5, Math.random()-0.5];
                let len  = Math.hypot(axis[0], axis[1], axis[2]);
                if (len === 0) { axis = [1,0,0]; len = 1; }
                axis = axis.map(a => a/len);
                const angle = (Math.random()-0.5) * scatterRad;
                const c = Math.cos(angle/2);
                const s = Math.sin(angle/2);
                quat[i*4+0] = axis[0]*s;
                quat[i*4+1] = axis[1]*s;
                quat[i*4+2] = axis[2]*s;
                quat[i*4+3] = c;

                const nose = rotateVecByQuat([0, 0, 1], [quat[i*4+0], quat[i*4+1], quat[i*4+2], quat[i*4+3]]);
                const speed = (Math.random() - 0.5 ) * speedSpread + speedBase;
                vel[i*4+0] = nose[0] * speed;
                vel[i*4+1] = nose[1] * speed;
                vel[i*4+2] = nose[2] * speed;
                vel[i*4+3] = 0.0;
            }
        }

        console.log("Init pos[0..7]:", Array.from(pos.slice(0, 8)));
        console.log("Init vel[0..7]:", Array.from(vel.slice(0, 8)));
        console.log("Init quat[0..7]:", Array.from(quat.slice(0, 8)));

        const texs = this.textures[this.pingPong];
        const texsNext = this.textures[1 - this.pingPong];
        
        const upload = (tex, data) => {
            gl.bindTexture(gl.TEXTURE_2D, tex);
            gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, this.size, this.size, gl.RGBA, gl.FLOAT, data);
        };

        upload(texs.pos, pos); upload(texsNext.pos, pos);
        upload(texs.quat, quat); upload(texsNext.quat, quat);
        upload(texs.vel, vel); upload(texsNext.vel, vel);
        upload(texs.angVel, zero); upload(texsNext.angVel, zero);
        upload(texs.force, zero); upload(texsNext.force, zero);
        upload(texs.torque, zero); upload(texsNext.torque, zero);
    }

    readback(verbosity) {
        if (verbosity <= 0) return;
        const gl = this.gl;
        const count = this.nParticles;
        const dataPos = new Float32Array(this.size * this.size * 4);
        const dataVel = new Float32Array(this.size * this.size * 4);
        
        gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbos[this.pingPong]);
        
        gl.readBuffer(gl.COLOR_ATTACHMENT0);
        gl.readPixels(0, 0, this.size, this.size, gl.RGBA, gl.FLOAT, dataPos);
        
        gl.readBuffer(gl.COLOR_ATTACHMENT1);
        gl.readPixels(0, 0, this.size, this.size, gl.RGBA, gl.FLOAT, dataVel);
        
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);

        let avgPos = [0,0,0];
        let avgVel = [0,0,0];
        
        for (let i = 0; i < count; i++) {
            avgPos[0] += dataPos[i*4+0];
            avgPos[1] += dataPos[i*4+1];
            avgPos[2] += dataPos[i*4+2];
            
            avgVel[0] += dataVel[i*4+0];
            avgVel[1] += dataVel[i*4+1];
            avgVel[2] += dataVel[i*4+2];
        }
        
        avgPos = avgPos.map(v => v / count);
        avgVel = avgVel.map(v => v / count);

        console.log(`[GPU Readback] Avg Pos: (${avgPos[0].toFixed(3)}, ${avgPos[1].toFixed(3)}, ${avgPos[2].toFixed(3)}) | Avg Vel: (${avgVel[0].toFixed(3)}, ${avgVel[1].toFixed(3)}, ${avgVel[2].toFixed(3)})`);
        
        if (verbosity >= 2) {
            console.log("First Particle Pos:", dataPos.slice(0, 3));
            console.log("First Particle Vel:", dataVel.slice(0, 3));
        }
    }

    step(dt, gravity, bodyDef) {
        if (!this.program) return;
        if (bodyDef) this.bodyDef = bodyDef;
        const gl = this.gl;
        const next = 1 - this.pingPong;
        const read = this.textures[this.pingPong];
        const write = this.fbos[next];

        gl.bindFramebuffer(gl.FRAMEBUFFER, write);
        gl.drawBuffers([gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1, gl.COLOR_ATTACHMENT2, gl.COLOR_ATTACHMENT3, gl.COLOR_ATTACHMENT4, gl.COLOR_ATTACHMENT5]);
        gl.viewport(0, 0, this.size, this.size);
        gl.useProgram(this.program);

        gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, read.pos); gl.uniform1i(this.locations.u_tex_pos, 0);
        gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, read.vel); gl.uniform1i(this.locations.u_tex_vel, 1);
        gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, read.quat); gl.uniform1i(this.locations.u_tex_quat, 2);
        gl.activeTexture(gl.TEXTURE3); gl.bindTexture(gl.TEXTURE_2D, read.angVel); gl.uniform1i(this.locations.u_tex_ang_vel, 3);
        gl.activeTexture(gl.TEXTURE4); gl.bindTexture(gl.TEXTURE_2D, read.force); gl.uniform1i(this.locations.u_tex_force, 4);
        gl.activeTexture(gl.TEXTURE5); gl.bindTexture(gl.TEXTURE_2D, read.torque); gl.uniform1i(this.locations.u_tex_torque, 5);

        gl.uniform1f(this.locations.u_dt, dt);
        gl.uniform4f(this.locations.u_gravity, 0, gravity, 0, 0);
        gl.uniform1i(this.locations.u_point_count, this.bodyDef.points.length / 3);
        gl.uniform3fv(this.locations.u_points, this.bodyDef.points);
        gl.uniform3f(this.locations.u_inertia_inv, ...this.bodyDef.inertiaInv);

        gl.bindVertexArray(this.quadVAO);
        gl.drawArrays(gl.TRIANGLES, 0, 6);
        
        // Cleanup state to avoid leaking draw buffers into Three.js
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        if (gl.drawBuffers) gl.drawBuffers([gl.BACK]);
        gl.bindVertexArray(null);
        gl.useProgram(null);
        this.pingPong = next;
    }

    getOutput() {
        return this.textures[this.pingPong];
    }

    debugReadback() {
        const gl = this.gl;
        const texs = this.textures[this.pingPong];
        const data = new Float32Array(4);
        
        gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbos[this.pingPong]);
        
        const check = (attachment, name) => {
            gl.readBuffer(attachment);
            gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, data);
            console.log(`[Debug 0,0] ${name}:`, data[0], data[1], data[2], data[3]);
        };

        check(gl.COLOR_ATTACHMENT0, "Pos ");
        check(gl.COLOR_ATTACHMENT1, "Vel ");
        check(gl.COLOR_ATTACHMENT2, "Quat");
        check(gl.COLOR_ATTACHMENT3, "AngV");
        check(gl.COLOR_ATTACHMENT4, "Forc");
        check(gl.COLOR_ATTACHMENT5, "Torq");

        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    }
}

