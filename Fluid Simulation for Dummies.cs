

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


void set_bnd(int b, float[] x, int N)
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




