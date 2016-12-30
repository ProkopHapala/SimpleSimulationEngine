
function makeTextSprite( message, parameters ){
	if ( parameters === undefined ) parameters = {};
	var  font  = parameters.hasOwnProperty("font")  ? parameters["font"]  : "Arial";
	var  size  = parameters.hasOwnProperty("size")  ? parameters["size"]  : 18;
	var color = parameters.hasOwnProperty("color") ? parameters["color"] : { r:0, g:0, b:0, a:1.0 };
	var canvas   = document.createElement('canvas');
	var context  = canvas.getContext('2d');
	context.font = "Bold " + size + "px " + font;
	var metrics = context.measureText( message );
	var textWidth = metrics.width;
	context.fillStyle = "rgba(" + color.r + "," + color.g + "," + color.b + "," + color.a + ")";
	context.fillText( message, 0, size);
	var texture = new THREE.Texture(canvas) 
	texture.needsUpdate = true;
	var spriteMaterial = new THREE.SpriteMaterial( { map: texture, useScreenCoordinates: false } );
	var sprite = new THREE.Sprite( spriteMaterial );
	sprite.scale.set(100,50,1.0);
	return sprite;	
}

function truss( n, sz, y0, w, matwire, mat ){
	//object = new THREE.Mesh( , matwire );
	//object.position.set( 0, 200, 0);
	//scene.add( object );
	var mesh;
	beam = new THREE.CylinderGeometry( w, w, 2*(n-1)*sz, 4, 1, false );
	mesh = new THREE.Mesh( beam, mat ); mesh.position.x = -sz; mesh.position.y = y0; scene.add( mesh );
	mesh = new THREE.Mesh( beam, mat ); mesh.position.x = +sz; mesh.position.y = y0; scene.add( mesh );
	mesh = new THREE.Mesh( beam, mat ); mesh.position.z = -sz; mesh.position.y = y0; scene.add( mesh );
	mesh = new THREE.Mesh( beam, mat ); mesh.position.z = +sz; mesh.position.y = y0; scene.add( mesh );
	block  = new THREE.OctahedronGeometry( sz, 0 );
	for ( var i = 0; i < n; i ++ ) {
	    mesh = new THREE.Mesh( block, matwire );
	    mesh.position.y = y0 + 2*(i-0.5*n+0.5)*sz;
	    //mesh.updateMatrix();
	    //mesh.matrixAutoUpdate = false;
	    scene.add( mesh );
    }
}

function trapezGeom( x0, y0, x1, y1, sym ){
    var geom = new THREE.Geometry(); 
    var y00=0.0,y10=0.0;
    if( sym ){ y00=-y0; y10=-y1; };
    geom.vertices.push(new THREE.Vector3(y00,x0,0));
    geom.vertices.push(new THREE.Vector3(y0 ,x0,0));
    geom.vertices.push(new THREE.Vector3(y10,x1,0));
    geom.vertices.push(new THREE.Vector3(y1 ,x1,0));
    geom.faces.push( new THREE.Face3( 0, 1, 2 ) );
    geom.faces.push( new THREE.Face3( 1, 2, 3 ) );
    
    geom.computeFaceNormals();
    //geom.computeVertexNormals();    // requires correct face normals
    //var object = new THREE.Mesh( geom, mat );
    //scene.addObject(object);
    return geom;
}

function BuildSpaceShip_1( scene ){

/*
	var geometry = new THREE.CylinderGeometry( 0, 10, 30, 16, 1 );
	var material = new THREE.MeshPhongMaterial( { color:0xffffff, shading: THREE.FlatShading } );
	for ( var i = 0; i < 500; i ++ ) {
		var mesh = new THREE.Mesh( geometry, material );
		mesh.position.x = ( Math.random() - 0.5 ) * 1000;
		mesh.position.y = ( Math.random() - 0.5 ) * 1000;
		mesh.position.z = ( Math.random() - 0.5 ) * 1000;
		mesh.updateMatrix();
		mesh.matrixAutoUpdate = false;
		scene.add( mesh );
	}
*/

//var material = new THREE.MeshLambertMaterial( { map: map, side: THREE.DoubleSide } );
	matwire      = new THREE.MeshBasicMaterial( { color: 0x808080, wireframe: true } );
	var material = new THREE.MeshLambertMaterial( { side: THREE.DoubleSide } );
	var material = new THREE.MeshPhongMaterial( { side: THREE.DoubleSide } );
	
	var geom;
	var spritey;
	
	var txColor = {r:0, g:200, b:0, a:1.0};

    //  barrel
    //                                                     
	//object = new THREE.Mesh( new THREE.CylinderGeometry( 5, 5, 500, 40, 5, {openEnded: true} ), material, );
	//object = new THREE.Mesh( new THREE.CylinderGeometry( {radiusTop: 5, radiusBottom: 5, height: 500, radiusSegments: 16, heightSegments: 5, openEnded: true} ), material );
	//object = new THREE.Mesh( new THREE.CylinderGeometry(radiusTop, radiusBottom, height, radiusSegments, heightSegments, openEnded, thetaStart, thetaLength);
	object = new THREE.Mesh( new THREE.CylinderGeometry(        4.0,          5.0, 800.0,              16,              1,      true ),  material );
	//object = new THREE.Mesh( new THREE.CylinderGeometry(      5.0,          5.0, 500.0,              16,              1,      false ),  material );
	object.position.set( 0, 0, 0 );
	scene.add( object );
	
	truss( 22, 16.0, -20, 1.0, matwire, material );
	
	
	// reactor
	//object = new THREE.Mesh( new THREE.CylinderGeometry( 25.0, 10.0, 100.0, 32, 1 ), material );
	geom = new THREE.CylinderGeometry( 30.0, 10.0, 100.0, 4, 1, false, Math.PI*0.25 ); 
	geom.shading = THREE.FlatShading;
	object = new THREE.Mesh( geom, material );
	object.position.set( 0, -400.0, 0 );
	object.scale.z = 0.5;
	scene.add( object );
	sprite = makeTextSprite( "reactor", { size: 36, color: txColor } ); sprite.position.set(50,-400,40); scene.add( sprite );
	
	// shield
	/*
	geom = new THREE.CylinderGeometry( 0.0, 60.0, 100.0, 4, 1, true );
	geom.shading = THREE.FlatShading;
	object = new THREE.Mesh( geom, material );
	object.scale.z = 0.25;
	*/
	//object = new THREE.Mesh( new THREE.PlaneGeometry( 80.0, 160.0 ), material );
	object = new THREE.Mesh( trapezGeom( -80, 40, 80, 48, true ), material );
	object.position.set( 0, +320.0, 0 );
	object.rotation.x = -0.15;
	scene.add( object );
	sprite = makeTextSprite( "Shield mirror", { size: 36, color: txColor } ); sprite.position.set(50,+320,40); scene.add( sprite );
	
	// Truss
	//object = new THREE.Mesh( new THREE.OctahedronGeometry( 75.0, 0 ), matwire );
	//object.position.set( 0, 200, 0);
	//scene.add( object );
	
    

    // radiator
	//object = new THREE.Mesh( new THREE.RingGeometry( 10.0, 200.0, 64, 1, 100.0, Math.PI * 2 ), material ); object.position.set( 0, -100, 0 ); object.rotation.x = 3.14159265359/2;
	object = new THREE.Mesh( trapezGeom( -350.0, 20.0, 300.0, 100.0, true ), material ); object.position.set( 0, 0, 0 );
	scene.add( object );
	sprite = makeTextSprite( "radiator", { size: 36, color: txColor } ); sprite.position.set(100,+100,40); scene.add( sprite );
	
	// Facilities
	object = new THREE.Mesh( new THREE.TorusGeometry( 32.0, 8.0, 16, 64 ), material ); object.position.set( 0, 200, 0 );
	//object.rotation.x = 3.14159265359/2;
	scene.add( object );
	sprite = makeTextSprite( "Facilities", { size: 36, color: txColor } ); sprite.position.set(50,+180,40); scene.add( sprite );
	
    // propelant tank
	object = new THREE.Mesh( new THREE.TorusGeometry( 40.0, 16.0, 16, 64 ), material ); object.position.set( 0, 0, 0 ); scene.add( object );
	object = new THREE.Mesh( new THREE.TorusGeometry( 30.0, 16.0, 16, 64 ), material ); object.position.set( 0, 0, 0 ); scene.add( object );
	sprite = makeTextSprite( "Propelant Tank", { size: 36, color: txColor } ); sprite.position.set(60,+0.0,40); scene.add( sprite );
	
	/*
	object = new THREE.Mesh( new THREE.TorusGeometry( 80.0, 10.0, 16, 64 ), material );
	object.position.set( 0, 0, 0 );
	object.rotation.y = 3.14159265359/2;
	scene.add( object );
	*/
	
	 //sprite = makeTextSprite( " Hello, ", { size: 24, color: txColor } ); sprite.position.set(-85,105,55); scene.add( sprite );
	
}
