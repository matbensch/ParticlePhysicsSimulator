import * as THREE from 'https://threejsfundamentals.org/threejs/resources/threejs/r127/build/three.module.js';
import {OrbitControls} from 'https://threejsfundamentals.org/threejs/resources/threejs/r127/examples/jsm/controls/OrbitControls.js';

//necessary variables for running three.js
let scene,camera,renderer,controls;

//universal constant values
const G = 6.67*Math.pow(10,-11);
const k = 8.99*Math.pow(10,9);
const mu = 4*Math.PI*Math.pow(10,-7);

//set extreme radius
const max_radius = 75;
const min_radius = 0.1;

//html elements
var gBox = document.getElementById("gravityBox");
var eBox = document.getElementById("electricBox");
var mBox = document.getElementById("magnetBox");


//the size of the 3D canvas
const range = 100;

//maximum velocity
const vmax = 100;

//returns the dot product of 3-d vector arrays a and b
function dot(a,b)
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//returns the cross product of 3-d vector arrays a and b
function cross(a,b)
{
    return [ a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]  ];
}

//returns the product of some scalar c and vector a
function mult(c,a)
{
    return [c*a[0],c*a[1],c*a[2]];
}

//returns the sum of vectors a and b
function add(a,b)
{
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2]];
}

//returns the difference of vectors a and b
function sub(a,b)
{
    return add(a,mult(-1,b));
}

//returns the magnitude of vector a
function mag(a)
{
    return Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

//class Particle is of type THREE.Mesh 
//it stores:
//  geometry (shape)
//  material (color)
//  radius r
//  mass m
//  charge q
//  position object (x,y,z)
//  x-velocity vx
//  y-velocity vy
//  z-velocity vz

class Particle extends THREE.Mesh
{
    constructor(geometry,material,r,m,q,x,y,z,vx,vy,vz)
    {
        super(geometry,material);
        this.isParticle = true;
        this.r = r
        this.m = m*1000000000000;
        this.q = q*100;
        this.position.set(x,y,z);
        this.vx=vx;
        this.vy=vy;
        this.vz=vz;
        this.ax = 0;
        this.ay = 0;
        this.az = 0;
    }

    //use the current acceleration, velocity, and position 
    //to calculate and set the new position of the particle after time change of dt
    move(dt)
    {
        this.vx = Math.max(-vmax,Math.min(vmax,this.vx+this.ax*dt));
        this.vy = Math.max(-vmax,Math.min(vmax,this.vy+this.ay*dt));
        this.vz = Math.max(-vmax,Math.min(vmax,this.vz+this.az*dt));

        let x = this.position.x + this.vx * dt;
        let y = this.position.y + this.vy * dt;
        let z = this.position.z + this.vz * dt;

        this.position.set(x,y,z);
    }

}

//prevents the particles from moving outside the box
function fixEdgeCollisions()
{
    for(let i = 0;i<scene.children.length;i++)
    {
        let p = scene.children[i];
        if(p.isParticle==null) continue;

        if(p.position.x+p.r>range) { p.vx *= -1; p.position.x=range-p.r }
        if(p.position.x-p.r<-range) { p.vx *= -1; p.position.x=-range+p.r }

        if(p.position.y+p.r>range) { p.vy *= -1; p.position.y=range-p.r }
        if(p.position.y-p.r<-range) { p.vy *= -1; p.position.y=-range+p.r }

        if(p.position.z+p.r>range) { p.vz *= -1; p.position.z=range-p.r }
        if(p.position.z-p.r<-range) { p.vz *= -1; p.position.z=-range+p.r }
    }
}

//prevents the collision of two particles
function fixParticleCollisions()
{
    for(let i = 0;i<scene.children.length;i++)
    {
        for(let j=i+1;j<scene.children.length;j++)
        {
            let p = scene.children[i];
            let q = scene.children[j];

            if(p.isParticle==null || q.isParticle==null) continue;

            let x1 = [p.position.x,p.position.y,p.position.z];
            let x2 = [q.position.x,q.position.y,q.position.z];
            let dx = sub(x1,x2 );

            if(mag(dx)<p.r+q.r)
            {
                let u1 = [p.vx,p.vy,p.vz];
                let u2 = [q.vx,q.vy,q.vz];
                
                let du = sub(u1,u2);
                let k = mult(1/mag(dx),dx);
                let a = (2 / (1/p.m+1/q.m) )*dot(k,du)

                let v1 = sub(u1,mult(a/p.m,k));
                let v2 = add(u2,mult(a/q.m,k));
                
                p.vx = v1[0];
                p.vy = v1[1];
                p.vz = v1[2];
                
                q.vx = v2[0];
                q.vy = v2[1];
                q.vz = v2[2];

                let overlap = p.r+q.r-mag(dx);
                
               let delPos1 = add(x1,mult( mag(u1)/(mag(u1)+mag(u2)) * overlap, k ) );
               let delPos2 = sub(x2,mult( mag(u2)/(mag(u1)+mag(u2)) * overlap, k ) );
                
                p.position.x = delPos1[0];
                p.position.y = delPos1[1];
                p.position.z = delPos1[2];

                q.position.x = delPos2[0];
                q.position.y = delPos2[1];
                q.position.z = delPos2[2];
                
                

                let qt = (p.q+q.q)/2;
                p.q = qt;
                q.q = qt;

                p.material = new THREE.MeshBasicMaterial({color:p.q<0?0x0000ff:0xff0000});
                q.material = new THREE.MeshBasicMaterial({color:q.q<0?0x0000ff:0xff0000});


            }
            
        }

    }
}

// applies the force of gravity between two particles and adds it to the acceleration
function applyGravity()
{
    for(let i=0;i<scene.children.length;i++)
    {
        for(let j=i+1;j<scene.children.length;j++)
        {
            let p = scene.children[i];
            let q = scene.children[j];

            if(p.isParticle==null || q.isParticle==null)
            {
                continue;
            }

            let dx = q.position.x-p.position.x;
            let dy = q.position.y-p.position.y;
            let dz = q.position.z-p.position.z;

            let dist = Math.max(0.01,Math.sqrt(dx*dx+dy*dy+dz*dz));

            let gravity = G * p.m * q.m / Math.pow(dist,2);


            let f_x = gravity * ( dx/dist );
            let f_y = gravity * ( dy/dist );
            let f_z = gravity * ( dz/dist );

            p.ax += f_x/p.m;
            p.ay += f_y/p.m;
            p.az += f_z/p.m;

            q.ax -= f_x/q.m;
            q.ay -= f_y/q.m;
            q.az -= f_z/q.m;
        }
    }
}

//applies the electrostatic force and adds to acceleration of particles
function applyElectricity()
{
    for(let i=0;i<scene.children.length;i++)
    {
        for(let j=i+1;j<scene.children.length;j++)
        {
            let p = scene.children[i];
            let q = scene.children[j];

            if(p.isParticle==null || q.isParticle==null)
            {
                continue;
            }

            let dx = q.position.x-p.position.x;
            let dy = q.position.y-p.position.y;
            let dz = q.position.z-p.position.z;

            let dist = Math.max(0.01,Math.sqrt(dx*dx+dy*dy+dz*dz));

            let electricity = k *Math.abs(p.q) * Math.abs(q.q) / Math.pow(dist,2);
            let attract = ( (p.q<0 && q.q<0)||(p.q>0 && q.q>0) )?-1:1;
            let f_x = attract * electricity * ( dx/dist );
            let f_y = attract * electricity * ( dy/dist );
            let f_z = attract * electricity * ( dz/dist );

            p.ax += f_x/p.m;
            p.ay += f_y/p.m;
            p.az += f_z/p.m;
            q.ax -= f_x/q.m;
            q.ay -= f_y/q.m;
            q.az -= f_z/q.m;
        }
    }
}

//applies the magnetic force and adds to the acceleration of all particles
function applyMagnetism()
{
    for(let i = 0;i<scene.children.length;i++)
    {
        for(let j = 0;j<scene.children.length;j++)
        {
            let p = scene.children[i];
            let q = scene.children[j];

            if(p.isParticle==null || q.isParticle==null) continue;
            if(i==j) continue;
            if(q.q==0) continue;


            let qvvec = [q.vx,q.vy,q.vz];
            let qrvec = [p.position.x-q.position.x,p.position.y-q.position.y,p.position.z-q.position.z];

            let pvvec = [p.vx,p.vy,p.vz];

            let B = mult(mu/(4*Math.PI)*q.q/Math.pow(mag(qrvec),3),cross(qvvec,qrvec));
            let F = mult(p.q, cross(pvvec,B) );

            if(F[0]==NaN) continue;
            if(F[1]==NaN) continue;
            if(F[2]==NaN) continue;
            

            p.ax += F[0]/p.m;
            p.ay += F[1]/p.m;
            p.az += F[2]/p.m;


        }
    }
}

//resets acceleration of all particles to zero
function refresh()
{
    for(let i=0;i<scene.children.length;i++)
    {
        scene.children[i].ax=0;
        scene.children[i].ay=0;
        scene.children[i].az=0;
        
    }
}

//adds new particle to the scene
function addParticle(scene,r,m,q,x,y,z,vx=0,vy=0,vz=0)
{
    let geometry = new THREE.SphereGeometry(r,32,32);
    let material = new THREE.MeshBasicMaterial({color: q>=0?0xff0000:0x0000ff});

    scene.add( new Particle(geometry,material,r,m,q,x,y,z,vx,vy,vz) );
}

//initializes all three.js values
function init()
{
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(75,window.innerWidth/(window.innerHeight),0.1,1000);
    
    renderer = new THREE.WebGLRenderer({antialias: true});
    renderer.setSize(window.innerWidth,window.innerHeight);
    renderer.setClearColor(0xffffff,1);

    controls = new OrbitControls(camera,renderer.domElement);
    camera.position.set(0,0,3*range);
    controls.update();

    var box = new THREE.BufferGeometry();
    var vertices = new Float32Array(
        [
            range,range,range,
            range,range,-range,
            range,-range,range,
            range,-range,-range,
            -range,range,range,
            -range,range,-range,
            -range,-range,range,
            -range,-range,-range
        ]
    )

    box.setAttribute('position',new THREE.BufferAttribute(vertices,3))

	//box.vertices.push( new THREE.Vector3( -range, -range, -range ) );
	//box.vertices.push( new THREE.Vector3( range, range, range ) );
	var boxMesh = new THREE.Line( box );
    scene.add( new THREE.BoxHelper( boxMesh, 'black' ) );
    
    var axes = new THREE.AxesHelper(100);
    axes.position.set(-100,-100,-100);
    scene.add(axes);

    var light = new THREE.DirectionalLight( 0x000000, .8 );
	light.position.set( -range, range, 0 );
	camera.add( light );
    scene.add( camera );
    
    var ambient = new THREE.AmbientLight( 0x555555 );
	scene.add( ambient );

    document.body.appendChild(renderer.domElement);
}

//changes size of scene with respect to window size
function onWindowResize()
{
    camera.aspect = window.innerWidth/(window.innerHeight);
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth,(window.innerHeight));

}
window.addEventListener('resize',onWindowResize,false);

//takes in user input of particle
function userParticle()
{
    let r = parseFloat(document.getElementById("radius").value )
    let m = parseFloat(document.getElementById("mass").value )
    let q = parseFloat(document.getElementById("charge").value )

    let x = parseFloat(document.getElementById("xpos").value )
    let y = parseFloat(document.getElementById("ypos").value )
    let z = parseFloat(document.getElementById("zpos").value )
    addParticle(scene,r,m,q,x,y,z);
}
document.getElementById("userParticleButton").onclick = userParticle;

//animates the scene
async function animate()
{
    requestAnimationFrame(animate);

    if(gBox.checked) applyGravity();
    if(eBox.checked) applyElectricity();
    if(mBox.checked) applyMagnetism();

    fixParticleCollisions();
    fixEdgeCollisions();

    for(let i = 0;i<scene.children.length;i++)
    {
        if(scene.children[i].isParticle==null) continue;
        scene.children[i].move(0.05);
    }

    refresh();

    controls.update();
    renderer.render(scene,camera);
}

init();
animate();