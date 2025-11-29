import { MeshRenderer } from '../../common_js/MeshRenderer.js';
import { Draw3D } from '../../common_js/Draw3D.js';
import { Logger } from '../../common_js/Logger.js';

export class MoleculeRenderer extends MeshRenderer {
    constructor(scene, system, shaders, mmParams) {
        super(scene, shaders, system.capacity);
        this.system = system;
        this.mmParams = mmParams;

        this.selectionMesh = null;
        this.axesHelper = null;

        this.labelMode = 'none';

        this.initSelection();
    }

    initSelection() {
        // --- Selection (InstancedMesh - Rings) ---
        this.selectionBaseGeo = Draw3D.createOctSphereGeometry(16, 1.0);
        this.selectionMesh = Draw3D.createTextureBasedSelectionLines(
            this.capacity,
            this.selectionBaseGeo,
            this.shaders.selection,
            {
                ...this.commonUniforms,
                uColor: { value: new THREE.Vector4(1.0, 1.0, 0.0, 0.8) } // Yellow
            }
        );
        this.scene.add(this.selectionMesh);
    }

    update() {
        // Main update loop
        if (this.system.isDirty) {
            this.updateStructure();
            this.updatePositions(); // Structure change implies position update too
            this.system.isDirty = false;
        }

        // Update Labels (visibility check handled in parent)
        // But we need to update content if mode changes or structure changes.
        // Here we just ensure visibility logic if needed.
    }

    updatePositions() {
        super.updatePositions(this.system.pos, this.system.nAtoms);
    }

    updateStructure() {
        const nAtoms = this.system.nAtoms;

        // Update Particles
        this.updateParticles(nAtoms, (i) => {
            const type = this.system.types[i];
            // Use MMParams
            let col = [1, 0, 1]; // Default magenta
            if (this.mmParams) {
                col = this.mmParams.getColor(type);
            }
            return col;
        }, (i) => {
            const type = this.system.types[i];
            let radius = 1.0;
            if (this.mmParams) {
                radius = this.mmParams.getRadius(type) * 0.4;
            }
            return radius;
        });

        // Update Bonds
        this.updateBonds(this.system.bonds);

        this.updateSelection();
        this.updateLabelsContent(); // Update labels when structure changes

        if (window.VERBOSITY_LEVEL >= Logger.INFO) {
            // window.logger.info(`Renderer structure updated: ${ nAtoms } atoms.`);
        }
    }

    updateSelection() {
        const selectedIDs = Array.from(this.system.selection);
        const count = selectedIDs.length;

        // Update instance count for InstancedBufferGeometry
        this.selectionMesh.geometry.instanceCount = count;

        const idAttr = this.selectionMesh.geometry.getAttribute('aAtomID');
        for (let i = 0; i < count; i++) {
            idAttr.setX(i, selectedIDs[i]);
        }
        idAttr.needsUpdate = true;
    }

    updateLabelsContent() {
        if (!this.labelMesh || this.labelMode === 'none') return;

        const nAtoms = this.system.nAtoms;

        this.updateLabels((i) => {
            if (this.labelMode === 'id') {
                return i.toString();
            } else if (this.labelMode === 'element') {
                const type = this.system.types[i];
                if (this.mmParams && this.mmParams.byAtomicNumber[type]) {
                    return this.mmParams.byAtomicNumber[type].name;
                }
                return type.toString();
            } else if (this.labelMode === 'type') {
                return "?"; // Placeholder as per previous logic
            }
            return "";
        }, nAtoms);
    }

    toggleAxes(visible) {
        if (!this.axesHelper) {
            this.axesHelper = new THREE.AxesHelper(5); // 5 units length
            this.scene.add(this.axesHelper);
        }
        this.axesHelper.visible = visible;
    }

    setLabelMode(mode) {
        if (this.labelMode !== mode) {
            this.labelMode = mode;
            if (this.labelMesh) {
                this.labelMesh.visible = (mode !== 'none');
                if (this.labelMesh.visible) {
                    this.updateLabelsContent();
                }
            }
        }
    }

    setLabelStyle(colorHex, scale) {
        if (this.labelMesh) {
            if (colorHex) {
                this.labelMesh.material.uniforms.uColor.value.set(colorHex);
            }
            if (scale !== undefined && scale !== null) {
                this.labelMesh.material.uniforms.uScale.value = parseFloat(scale);
            }
        }
    }
}

