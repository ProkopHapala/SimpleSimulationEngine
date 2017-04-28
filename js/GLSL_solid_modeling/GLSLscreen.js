	// ===================== 
    //   THREE JS MAIN
    // =====================
    // global variables
    var renderer;
    var scene;
    var camera;

    var uniforms = {};
    var control;
    var mousePosOld;

	var basicShader;
	var mesh;
    
    var camQuat = new  THREE.Quaternion();
    //var camMat  = new THREE.Matrix3();
    
    function handleMouseMove(event) {
        var dot, eventDoc, doc, body, pageX, pageY;        
        if (mousePosOld) {
            var dx = (event.clientX-mousePosOld.x)*1.0;
            var dy = (event.clientY-mousePosOld.y)*1.0;
            //console.log( "x: "+ dx +" y: "+ dy );
            //var v = new THREE.Vector3( dx, dy, 0.0 );//.normalize();
            //var q = new THREE.Quaternion().setFromEuler( v );
            //THREE.Quaternion.slerp( camQuat, q, camQuat, 0.07 );
            //var q = new THREE.Quaternion( 0.01, 0.01, 0.01, 1.0 );
            var q = new THREE.Quaternion( -dy*0.002, 0.0, dx*0.002, 1.0 );
            camQuat.multiply(q).normalize();
            //console.log( q.x+" "+ q.y+" "+ q.z+" "+ q.w +" "+ camQuat.x+" "+ camQuat.y+" "+ camQuat.z+" "+ camQuat.w  );
            
        }else{
            mousePosOld = new THREE.Vector3();
        }
        mousePosOld.x = event.clientX;
        mousePosOld.y = event.clientY;
        
    }
    
    function render() {
        //console.log("x: " + event.clientX + ", y: " + event.clientY );
        
        var camMat4 = new THREE.Matrix4(); camMat4.compose ( new THREE.Vector3(0.0,0.0,0.0), camQuat, new THREE.Vector3(1.0,1.0,1.0) );
        var camMat_ = new THREE.Matrix3(); camMat_.getNormalMatrix ( camMat4 );
        //console.log( camMat_.elements +" "+ camQuat.x+" "+ camQuat.y+" "+ camQuat.z+" "+ camQuat.w  );
        uniforms.camMat.value = camMat_;
        
        renderer.render(scene, camera);
        uniforms.time.value += 0.05;
        requestAnimationFrame(render);
    }
    
    function updateShader(element){
		console.log("updateShader");
		console.log(element.value);

		var shader_code="";
		//jQuery.get('Primitives.glsl', function(data){ shader_code += data; } );
		//console.log(document.getElementById("txtRayTracer"));
		// see:
		// http://stackoverflow.com/questions/6348207/making-a-paragraph-in-html-contain-a-text-from-a-file
		// http://stackoverflow.com/questions/36659202/read-data-in-html-object-tag
		
		
		shader_code += document.getElementById("txtPrimitives").contentDocument.body.childNodes[0].textContent;
		shader_code += element.value;
		shader_code += document.getElementById("txtRayTracer").contentDocument.body.childNodes[0].textContent;
		
		//shader_code += document.getElementById("txtPrimitives").contentDocument.getElementsByTagName('body')[0].innerHTML;
		
		
		
		//shader_code = document.getElementById('Primitives').text + element.value + document.getElementById('RayTracer').text;
		//console.log("============");
		//console.log("============");
		//console.log("============");
		//console.log(shader_code);
	    //console.log("============");
		//console.log("============");
		//console.log("============");
		
        var material = new THREE.ShaderMaterial({
            uniforms: uniforms,
            vertexShader: basicShader.vertexShader,
            fragmentShader: shader_code,
        });
        		
		try {
		    var mesh = scene.getObjectByName( "Mesh1" );
            scene.remove( mesh );
        }catch(err) {
            console.log( " updateShader cannot remove Mesh1" );
        }; 

		mesh = new THREE.Mesh( new THREE.PlaneBufferGeometry( 200, 200 ), material );
		mesh.name = "Mesh1";
		scene.add( mesh );
	}
	
	function init_GLSLScreen(){
	
	}

    function init_GLSLScreen(){
        document.onmousemove = handleMouseMove;
    
        scene    = new THREE.Scene();
        camera   = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
        renderer = new THREE.WebGLRenderer();
        renderer.setClearColor(0x000000, 1.0);
        renderer.setSize(window.innerWidth, window.innerHeight);

        // SHADER
        basicShader = THREE.ShaderLib['normal'];
		uniforms = {
				time      : { value: 1.0 },
				resolution: { value: new THREE.Vector2() },
				camMat    : { value: new THREE.Matrix3() }
		};
		uniforms.resolution.value.x = renderer.domElement.width;
		uniforms.resolution.value.y = renderer.domElement.height;
        updateShader( document.getElementById("txtScene") );
				
        camera.position.x = 0.0;
        camera.position.y = 0.0;
        camera.position.z = 100.0;
        camera.lookAt(scene.position);
        
        //document.body.appendChild(renderer.domElement);
        document.getElementById("divLeft").appendChild( renderer.domElement );
        control = new function () {
            this.rotationSpeed = 0.005;
            this.scale = 1;
        };
        //addControls(control);
        // call the render function
        render();
    }

    // calls the init function when the window is done loading.
    //window.onload = init;
