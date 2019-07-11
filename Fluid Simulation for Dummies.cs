using System;

int IX(int x, int y, int z,int N)
{
	return x+y*N+z*N*N
}

struct FluidCube
{
	int size;
	float dt;
	float diff;
	float visc;
	
	float[] s;
	float[] density;
	
	float[] Vx;
	float[] Vy;
	float[] Vz;
	
	float[] Vx0;
	float[] Vy0;
	float[] Vz0;
}

// add density to a cell
void FluidCubeAddDensity(FluidCube cube, int x, int y, int z, float amout)
{
	int N = cube.size;
	cube.density[IX(x,y,z,N)] += amount;
}

// add velocity to a cell
void FluidCubeAddVelocity(FluidCube cube, int x, int y, int z, float amountX, float amountY, float amountZ)
{
	int N = cube.size;
	int index = IX(x,y,z,N);
	
	cube.Vx[index] += amountX;
	cube.Vy[index] += amountY;
	cube.Vz[index] += amountZ;
}

// set up the fluid cube object
FluidCube FluidCubeCreate(int size, int diffusion, int viscosity, float dt)
{
	FluidCube cube;
	int N = size;
	
	cube.size = size;
	cube.dt = dt;
	cube.diff = diffusion;
	cube.visc = viscosity;
	
	cube.s = new float[N*N*N];
	cube.density = new float[N*N*N];
	
	cube.Vx = new float[N*N*N];
	cube.Vy = new float[N*N*N];
	cube.Vz = new float[N*N*N];
	
	cube.Vx0 = new float[N*N*N];
	cube.Vy0 = new float[N*N*N];
	cube.Vz0 = new float[N*N*N];
	
	return cube;
}


float[] set_bnd(int b, float[] x, int N)
{
	for(int j = 1; j < N - 1; j++) 
	{
        for(int i = 1; i < N - 1; i++) 
		{
			if(b == 3)
			{
				x[IX(i, j, 0, N)] = -x[IX(i, j, 1,N)];
			}
			else
			{
				x[IX(i, j, 0, N)] = x[IX(i, j, 1,N)];
			}
			
            if(b == 3)
			{
				x[IX(i, j, N-1, N)] = -x[IX(i, j, N-2, N)];
			}
			else
			{
				x[IX(i, j, N-1, N)] = x[IX(i, j, N-2, N)];
			}
        }
    }
	
	for(int k = 1; k < N - 1; k++)
	{
        for(int i = 1; i < N - 1; i++)
		{
			if(	b == 2)
			{
				x[IX(i, 0  , k, N)] = -x[IX(i, 1  , k, N)];
			}
			else
			{
				x[IX(i, 0  , k, N)] = x[IX(i, 1  , k, N)];
			}
			if(	b == 2)
			{
				x[IX(i, N-1, k, N)] = -x[IX(i, N-2, k, N)];
			{
			else
			{		
				x[IX(i, N-1, k, N)] = x[IX(i, N-2, k, N)];
			}
        }
    }
	
	for(int k = 1; k < N - 1; k++)		
	{
        for(int j = 1; j < N - 1; j++) 
		{
			if(	b == 1)
			{
				x[IX(0  , j, k, N)] = -x[IX(1  , j, k, N)];
			}
			else
			{
				x[IX(0  , j, k, N)] = x[IX(1  , j, k, N)];
			}
			if(	b == 1)
			{
				x[IX(N-1, j, k, N)] = -x[IX(N-2, j, k, N)];
			}
			else
			{
				x[IX(N-1, j, k, N)] =  x[IX(N-2, j, k, N)];
			}
        }
    }
	
	x[IX(0, 0, 0, N)]       = 0.33f * (x[IX(1, 0, 0, N)]
                                  + x[IX(0, 1, 0, N)]
                                  + x[IX(0, 0, 1, N)]);
    x[IX(0, N-1, 0, N)]     = 0.33f * (x[IX(1, N-1, 0, N)]
                                  + x[IX(0, N-2, 0, N)]
                                  + x[IX(0, N-1, 1, N)]);
    x[IX(0, 0, N-1, N)]     = 0.33f * (x[IX(1, 0, N-1, N)]
                                  + x[IX(0, 1, N-1, N)]
                                  + x[IX(0, 0, N-2, N)]);
    x[IX(0, N-1, N-1, N)]   = 0.33f * (x[IX(1, N-1, N-1, N)]
                                  + x[IX(0, N-2, N-1, N)]
                                  + x[IX(0, N-1, N-2, N)]);
    x[IX(N-1, 0, 0, N)]     = 0.33f * (x[IX(N-2, 0, 0, N)]
                                  + x[IX(N-1, 1, 0, N)]
                                  + x[IX(N-1, 0, 1, N)]);
    x[IX(N-1, N-1, 0, N)]   = 0.33f * (x[IX(N-2, N-1, 0, N)]
                                  + x[IX(N-1, N-2, 0, N)]
                                  + x[IX(N-1, N-1, 1, N)]);
    x[IX(N-1, 0, N-1, N)]   = 0.33f * (x[IX(N-2, 0, N-1, N)]
                                  + x[IX(N-1, 1, N-1, N)]
                                  + x[IX(N-1, 0, N-2, N)]);
    x[IX(N-1, N-1, N-1, N)] = 0.33f * (x[IX(N-2, N-1, N-1, N)]
                                  + x[IX(N-1, N-2, N-1, N)]
                                  + x[IX(N-1, N-1, N-2, N)]);
								  
	return x;
}

float[] lin_solve(int b, float[] x, float[] x0, float a, float c, int iter, int N)
{
	float cRecip = 1.0f / c;
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m, N)] =
                        (x0[IX(i, j, m, N)]
                            + a*(    x[IX(i+1, j  , m  , N)]
                                    +x[IX(i-1, j  , m  , N)]
                                    +x[IX(i  , j+1, m  , N)]
                                    +x[IX(i  , j-1, m  , N)]
                                    +x[IX(i  , j  , m+1, N)]
                                    +x[IX(i  , j  , m-1, N)]
                           )) * cRecip;
                }
            }
        }
		x = set_bnd(b, x, N);
	}
	return x;
}

float[] diffuse(int b, float[] x, float[] x0, float diff, float dt, int iter, int N)
{
	float a = dt * diff * (N - 2) * (N - 2);
	x = lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
	return x;
}

float[] advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float[] velocZ, float dt, int N)
{
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    float dtz = dt * (N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y, z;
    
    float Nfloat = N;
    float ifloat, jfloat, kfloat;
    int i, j, k;
	
	
    for(k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j, k, N)];
                tmp2 = dty * velocY[IX(i, j, k, N)];
                tmp3 = dtz * velocZ[IX(i, j, k, N)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                z    = kfloat - tmp3;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = Math.Floor(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = Math.Floor(y);
                j1 = j0 + 1.0f; 
                if(z < 0.5f) z = 0.5f;
                if(z > Nfloat + 0.5f) z = Nfloat + 0.5f;
                k0 = Math.Floor(z);
                k1 = k0 + 1.0f;
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                u1 = z - k0;
                u0 = 1.0f - u1;
                
                int i0i = i0;
                int i1i = i1;
                int j0i = j0;
                int j1i = j1;
                int k0i = k0;
                int k1i = k1;
                
                d[IX(i, j, k, N)] = 
                
                    s0 * ( t0 * (u0 * d0[IX(i0i, j0i, k0i, N)]
                                +u1 * d0[IX(i0i, j0i, k1i, N)])
                        +( t1 * (u0 * d0[IX(i0i, j1i, k0i, N)]
                                +u1 * d0[IX(i0i, j1i, k1i, N)])))
                   +s1 * ( t0 * (u0 * d0[IX(i1i, j0i, k0i, N)]
                                +u1 * d0[IX(i1i, j0i, k1i, N)])
                        +( t1 * (u0 * d0[IX(i1i, j1i, k0i, N)]
                                +u1 * d0[IX(i1i, j1i, k1i, N)])));
            }
        }
    }
	d = set_bnd(b, d, N);
	return d;
}

void projet(ref float[] velocX, ref float[] velocY, ref float[] velocZ, float[] p, float[] div, int iter, int N)
{
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j, k, N)] = -0.5f*(
                         velocX[IX(i+1, j  , k  , N)]
                        -velocX[IX(i-1, j  , k  , N)]
                        +velocY[IX(i  , j+1, k  , N)]
                        -velocY[IX(i  , j-1, k  , N)]
                        +velocZ[IX(i  , j  , k+1, N)]
                        -velocZ[IX(i  , j  , k-1, N)]
                    )/N;
                p[IX(i, j, k, N)] = 0;
            }
        }
    }
    div = set_bnd(0, div, N); 
    p = set_bnd(0, p, N);
    p = lin_solve(0, p, div, 1, 6, iter, N);
    
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j, k, N)] -= 0.5f * (  p[IX(i+1, j, k, N)]
                                                -p[IX(i-1, j, k, N)]) * N;
                velocY[IX(i, j, k, N)] -= 0.5f * (  p[IX(i, j+1, k, N)]
                                                -p[IX(i, j-1, k, N)]) * N;
                velocZ[IX(i, j, k, N)] -= 0.5f * (  p[IX(i, j, k+1, N)]
                                                -p[IX(i, j, k-1, N)]) * N;
            }
        }
    }
   velocX = set_bnd(1, velocX, N);
   velocY = set_bnd(2, velocY, N);
   velocZ = set_bnd(3, velocZ, N);
	
}

FluidCube FluidCubeStep(FluidCube cube)
{
	/*
    int N          = cube.size;	
    float visc     = cube.visc;	
    float diff     = cube.diff;
    float dt       = cube.dt; 	
    float Vx      = cube.Vx;
    float Vy      = cube.Vy;
    float Vz      = cube.Vz;
    float Vx0     = cube.Vx0;
    float Vy0     = cube.Vy0;
    float Vz0     = cube.Vz0;
    float s       = cube.s;
    float density = cube.density;
	*/
    
    cube.Vx0 = diffuse(1, cube.Vx0, cube.Vx, cube.visc, cube.dt, 4, cube.size);
    cube.Vy0 = diffuse(2, cube.Vy0, cube.Vy, cube.visc, cube.dt, 4, cube.size);
    cube.Vz0 = diffuse(3, cube.Vz0, cube.Vz, cube.visc, cube.dt, 4, cube.size);
    
    project(ref cube.Vx0, ref cube.Vy0, ref cube.Vz0, cube.Vx, cube.Vy, 4, cube.size);
    
    cube.Vx = advect(1, cube.Vx, cube.Vx0, cube.Vx0, cube.Vy0, cube.Vz0, cube.dt, cube.size);
    cube.Vy = advect(2, cube.Vy, cube.Vy0, cube.Vx0, cube.Vy0, cube.Vz0, cube.dt, cube.size);
    cube.Vz = advect(3, cube.Vz, cube.Vz0,cube.Vx0, cube.Vy0, cube.Vz0, cube.dt, cube.size);
    
    project(ref cube.Vx, ref cube.Vy, ref cube.Vz, cube.Vx0, cube.Vy0, 4, cube.size);
    
    cube.s = diffuse(0, cube.s, cube.density, cube.diff, cube.dt, 4, cube.size);
    cube.density = advect(0, cube.density, cube.s, cube.Vx, cube.Vy, cube.Vz, cube.dt, cube.size);
	
	return cube;
}




