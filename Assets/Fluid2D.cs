    using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading;
using UnityEngine;
using Debug = UnityEngine.Debug;

namespace DefaultNamespace
{
    public class Fluid2D
    {
        public static double MAX_SPEED = 9.2;
        private double height;
        private int numX;
        private int numY;
        double[,] u;
        double[,] v;
        public int[,] s;
        public double[,] smokeField;
        private double overrelaxation = 1.99;
        private string mode;
        public double[] speed;
        static double[,] auxU;
        static double[,] auxV;
        static double[,] auxW;
        public double[,] p;
        private int separation;
        private ComputeShader computeShader;
        public double density;
        private int numberOfFans;

        public Fluid2D(double density ,int numX, int numY, float height, double[] speed, string mode,ComputeShader computeShader,int separation,int numberOfFans)
        {
            this.density = density;
            this.numberOfFans = numberOfFans;
            this.computeShader = computeShader;
            this.mode = mode;
            this.speed = speed;
            /*for (int i = 0; i < numberOfFans; i++)
            {
                this.speed[i] = speed[i];
            }
            Array.Reverse(this.speed);*/
            
            this.height = height;
            this.separation = separation;
            this.numX = numX + 2; //para añadir las celdas de los bordes
            this.numY = numY + 2;
            this.u = new double[this.numX, this.numY];
            this.v = new double[this.numX, this.numY];
            this.p = new double[this.numX, this.numY];
            this.s = new int[this.numX, this.numY];
            this.smokeField = new double[this.numX, this.numY];
            for (int i = 0; i < this.numX; i++) //we want to have limits, so we put 0 in the borders an 1 inside the field
            {
                for (int j = 0; j < this.numY; j++)
                {
                    
                        this.s[i, j] = 0;
                        this.smokeField[i, j] = MAX_SPEED;
                        this.u[i, j] = 0;
                        this.v[i, j] = 0;
                        this.p[i, j] = 0;


                }
            }

            for (int i = 1;
                 i < this.numX - 1;
                 i++) //we want to have limits, so we put 0 in the borders an 1 inside the field
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    this.s[i, j] = 1;
                    
                }
            }


            //this.w = new float[this.numX,this.numY,this.numZ];
        }

        public void setSpeed(double[] speed)
        {
            this.speed = new double[numberOfFans];
            for (int i = 0; i < numberOfFans; i++)
            {
                this.speed[i] = speed[i];
                
            }
            Array.Reverse(this.speed);
        }

        public void simulate(double dt)
        {
            
            if (mode.Equals("y0"))
            {
                for (int i = numX/2;
                     i < numX/2+2;
                     i++) //we want to have limits, so we put 0 in the borders an 1 inside the field
                {
                    
                        v[i, 1] = speed[i];
                        //smokeField[i, 1, k] = 0;
                   
                }
                for (int i = 0; i < numX; i++)
                {
                    this.s[i, numY - 1] = 1;
                    
                }
            }
            else if (mode.Equals("y1"))
            {
            }
            else if (mode.Equals("x0"))
            {
                
                for (int j = 1;
                     j < numberOfFans+1;
                     j++) //we want to have limits, so we put 0 in the borders an 1 inside the field
                {
                    
                    u[1 ,j*separation-1] = speed[numberOfFans-(j)];
                    
                    //smokeField[i, 1, k] = 0;
                   
                }
                for (int i = 2; i < numY-1; i++)
                {
                    s[numX-1,i] = 1;
                }

                /*this.s[numX / 2, numY / 2] = 0;
                this.s[numX / 2, numY / 2-1] = 0;*/

                //s[numX / 2, numY / 2] = 0;
            }
            else if (mode.Equals("x1"))
            {
            }
            modifyVelocity(dt, -3);
            for (int i = 0; i < numX;i++)
            {
                for (int j = 0;j < numY;j++)
                {
                    p[i, j] = 0;
                }

            }
            forceIncomprensibility(200, dt);
            advection(dt);
            smoke(dt);
        }
        
         public void simulateGPU(double dt)
        {
            if (mode.Equals("y0"))
            {
                for (int i = numX/2;
                     i < numX/2+2;
                     i++) //we want to have limits, so we put 0 in the borders an 1 inside the field
                {

                    if (s[i, 1] != 0)
                    {
                        v[i, 1] = speed[i];
                        //smokeField[i, 1, k] = 0;
                    }
                }
                for (int i = 0; i < numX; i++)
                {
                    this.s[i, numY - 1] = 1;
                    
                }
            }
            else if (mode.Equals("x0"))
            {
                for (int j = 1;
                     j < numberOfFans+1;
                     j++) //we want to have limits, so we put 0 in the borders an 1 inside the field
                {
                    
                    u[1 ,j*separation-1] = speed[j-1];
                    
                    //smokeField[i, 1, k] = 0;
                   
                }
                for (int i = 2; i < numY-1; i++)
                {
                    s[numX-1,i] = 1;
                }

                //s[numX / 2, numY / 2] = 0;
            }
            else if (mode.Equals("x1"))
            {
            }
            modifyVelocity(dt, -3);
            for (int i = 0; i < numX;i++)
            {
                for (int j = 0;j < numY;j++)
                {
                    p[i, j] = 0;
                }

            }
            forceIncomprensibilityGPU(200, dt);
            advection(dt);
            smoke(dt);
        }
        

        public void modifyVelocity(double dt, double aceleration)
        {
           
         
            //first we need to take into account the base values for the different vectors, adding the gravity to the v vector
            for (int i = 1; i < this.numX; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    if (this.s[i, j] != 0.0 && this.s[i, j - 1] != 0.0)
                        this.v[i, j] += aceleration * 0.02;
                }
            }
      
        }
        
        public void forceIncomprensibilityGPU(int iterations, double dt)
        {
            double[] auxU = new double[this.numX * this.numY],
                auxV = new double[this.numX * this.numY],
                ubuffer = new double[this.numX * this.numY],
                vbuffer = new double[this.numX * this.numY],
                pbuffer = new double[this.numX * this.numY];
            int[] sbuffer = new int[this.numX * this.numY];
     
         
            int l;
            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    ubuffer[i  * numY + j  ] = this.u[i,j];
                    vbuffer[i  * numY + j  ] = this.v[i,j];
                    sbuffer[i  * numY + j  ] = this.s[i,j];
                    
                }
            }
            
            computeShader.SetInt("numY",this.numY);
            computeShader.SetInt("numX",this.numX);
            computeShader.SetFloat("dt",(float)dt);
            computeShader.SetFloat("density",(float)this.density);
            computeShader.SetFloat("height",(float)this.height);
            
            ComputeBuffer uComputeBuffer = new ComputeBuffer(ubuffer.Length,sizeof(double));
            ComputeBuffer vComputeBuffer = new ComputeBuffer(vbuffer.Length,sizeof(double));
            ComputeBuffer pComputeBuffer = new ComputeBuffer(pbuffer.Length,sizeof(double));
            ComputeBuffer auxUComputeBuffer = new ComputeBuffer(auxU.Length,sizeof(double));
            ComputeBuffer auxVComputeBuffer = new ComputeBuffer(auxV.Length,sizeof(double));
            ComputeBuffer sComputeBuffer = new ComputeBuffer(sbuffer.Length,sizeof(int));
            auxUComputeBuffer.SetData(ubuffer);
            auxVComputeBuffer.SetData(vbuffer);
            
            sComputeBuffer.SetData(sbuffer);
            
            computeShader.SetBuffer(0,"auxU",auxUComputeBuffer);
            computeShader.SetBuffer(0,"auxV",auxVComputeBuffer);
            computeShader.SetBuffer(0,"s",sComputeBuffer);
            
            
            
            for (l = 0; (l < iterations); l++)
            {
                uComputeBuffer.SetData(ubuffer);// fill the buffer with the new data
                vComputeBuffer.SetData(vbuffer);// 
                pComputeBuffer.SetData(pbuffer);
                computeShader.SetBuffer(0,"u",uComputeBuffer);
                computeShader.SetBuffer(0,"v",vComputeBuffer);
                computeShader.SetBuffer(0,"p",pComputeBuffer);
                computeShader.Dispatch(0,numX/10,numY/10,1);
                auxUComputeBuffer.GetData(ubuffer);//u = auxU
                auxVComputeBuffer.GetData(vbuffer);//v = auxV
                pComputeBuffer.GetData(pbuffer);
                
                
            }
         
            uComputeBuffer.Dispose();
            vComputeBuffer.Dispose();
            pComputeBuffer.Dispose();
            auxUComputeBuffer.Dispose();
            auxVComputeBuffer.Dispose();

            sComputeBuffer.Dispose();

            for (int i = 0; i < numX; i++)
            {
                for (int j = 0; j < numY; j++)
                {
                    
                        this.u[i, j] = ubuffer[i*numY+j];
                        this.v[i,j] =vbuffer[i*numY+j];
                        this.p[i,j] = pbuffer[i*numY+j];
                    
                }
            }
            String results = "";
            for (int i = 0; i < numX;i++)
            {
                for (int j = 0;j < numY;j++)
                {
                    results += " " + ubuffer[i * numY + j ];
                    
                }
                results += "\n";
                
            }
            
            /*Debug.Log(results);*/
            Debug.Log("iterations:"+ l);

        }

       
        public void forceIncomprensibility(int iterations, double dt)
        {
            double outflow = 0;
            int actualS = 0;
            double[,] outflowMat = new double[this.numX, this.numY];
            bool converged = false;
            double error = 0.00001;
            int l = 0;
            
            Stopwatch sw = new Stopwatch(); 
            sw.Start(); 



            
            double outflowAvg = 0;
            for (l = 0; (l < iterations && !converged); l++)
            {
                
                converged = true;
                outflowAvg = 0;
                for (int i = 1; i < this.numX - 1; i++)
                {
                    for (int j = 1; j < this.numY - 1; j++)
                    {

                        if (this.s[i, j] != 0)
                        {
                            actualS = this.s[i + 1, j] + this.s[i - 1, j] + this.s[i, j + 1] +
                                      this.s[i, j - 1];
                            if (actualS != 0)
                            {
                                outflow = this.u[i + 1, j] - this.u[i, j] + this.v[i, j + 1] -
                                          this.v[i, j];
                                outflow *= overrelaxation;
                                
                                u[i, j] = this.u[i, j] + outflow * this.s[i - 1, j] / actualS;
                                u[i + 1, j] =
                                    this.u[i + 1, j] - outflow * this.s[i + 1, j] / actualS;
                                
                                v[i, j] = this.v[i, j] + outflow * this.s[i, j - 1] / actualS;
                                v[i, j + 1] =
                                    this.v[i, j + 1] - outflow * this.s[i, j + 1] / actualS;

                               
                                converged = converged && Math.Abs(outflowMat[i, j] - outflow) < error;
                                outflowMat[i, j] = outflow;
                                outflowAvg += outflow;
                                this.p[i, j] += (-outflow / actualS) * (density * height / dt);
                            }
                        }

                    }
                }

            }
            //Debug.Log("iterations:"+ l); 
            /*String results = "";
            for (int i = 0; i < numX;i++)
            {
                for (int j = 0;j < numY;j++)
                {
                    results += " " + u[i , j ];
                    
                }
                results += "\n";
                
            }*/
            /*Debug.Log(results);*/
            sw.Stop(); 
            //Debug.Log("Time elapsed: {0}"+ sw.Elapsed.ToString("hh\\:mm\\:ss\\.fff"));
           
            for (int j = 0; j < this.numY; j++)
            {
                this.v[0, j] = this.v[1, j];
                this.v[numX-1, j] = this.v[numX-2, j];
            }
            for (int i = 0; i < this.numX; i++)
            {
                this.u[i, 0] = this.u[i, 1];
                this.u[i, numY-1] = this.u[i, numY-2];
            }
            
        }

        public void advection(double dt)
        {
            double x, y,  x0, y0;
            double u, v;
            int i0, j0, i1, j1;
            double[,] auxu = new double[this.numX, this.numY];
            double[,] auxv = new double[this.numX, this.numY];
            double b000, b001, b010, b011;
            for (int i = 1; i < this.numX - 1; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    
                        //u
                        if (this.s[i, j] != 0.0)
                        {
                            u = this.u[i, j];
                            v = (this.v[i, j] + this.v[i - 1, j] +
                                 this.v[i, j + 1] + this.v[i - 1, j + 1]) / 4;
                            
                            x = i * this.height - dt * u;
                            y = j * this.height + 0.5 * this.height - dt * v;
                            x = Math.Max(Math.Min(x, this.numX * height),
                                height); //we have to be sure to not go out of the field
                            y = Math.Max(Math.Min(y, this.numY * height), height);
               

                            x0 = 0;
                            y0 = this.height * 0.5;

                            i0 = (int)Math.Min(Math.Floor(x / height), this.numX - 1); //we calculate the position relative to the vectors but in cells to know the cell from wich we must get the velocity but taking care not to go out of the field
                            j0 = (int)Math.Min(Math.Floor((y - y0) / height), this.numY - 1);
                           

                            i1 = Math.Min(i0 + 1, this.numX - 1); //be sure the next cell is not out of the field either
                            j1 = Math.Min(j0 + 1, this.numY - 1);

                            b001 = (x - i0 * height) /
                                   height; //we calculate the position of each component inside the cell to make the weights
                            b011 = ((y - y0 )- j0 * height) / height;
                            b000 = 1 - b001;
                            b010 = 1 - b011;
                            auxu[i, j] = this.u[i0, j0] * b000 * b010 +
                                         this.u[i1, j0] * b001 * b010 +
                                         this.u[i0, j1] * b000 * b011 +
                                         this.u[i1, j1] * b001 * b011;

                        }


                        //v
                        if (this.s[i, j] != 0.0)
                        {
                            v = this.v[i, j];
                            u = (this.u[i, j - 1] + this.u[i, j] +
                                 this.u[i + 1, j - 1] + this.u[i + 1, j]) / 4;
                            

                            x = i * this.height + 0.5 * this.height - dt * u;
                            y = j * this.height - dt * v;
                            x = Math.Max(Math.Min(x, this.numX * height),
                                height); //we have to be sure to not go out of the field
                            y = Math.Max(Math.Min(y, this.numY * height), height);
                           

                            x0 = this.height * 0.5;
                            y0 = 0;

                            i0 = (int)Math.Min(Math.Floor((x - x0) / height),
                                this.numX -
                                1); //we calculate the position relative to the vectors but in cells to know the cell from wich we must get the velocity but taking care not to go out of the field
                            j0 = (int)Math.Min(Math.Floor(y / height), this.numY - 1);
                          

                            i1 = Math.Min(i0 + 1,
                                this.numX - 1); //be sure the next cell is not out of the field either
                            j1 = Math.Min(j0 + 1, this.numY - 1);
                     
                            b001 = (x - x0 - i0 * height) /
                                   height; //we calculate the position of each component inside the cell to make the weights
                            b011 = (y - j0 * height) / height;
                            
                            b000 = 1 - b001;
                            b010 = 1 - b011;


                            auxv[i, j] = this.v[i0, j0] * b000 * b010 +
                                         this.v[i1, j0] * b001 * b010 
                                         +this.v[i0, j1] * b000 * b011 +
                                         this.v[i1, j1] * b001 * b011;

                        }

                    


                }
            }

            this.v = auxv;
            this.u = auxu;

        }

        public void smoke(double dt)
        {
            double u, v;
            double x, y;
            int cell0x;
            int cell0y;
            double[,] auxSmokeField = new double[this.numX, this.numY];
            
           
            for (int i = 1; i < this.numX-1; i++)
            {
                for (int j = 1; j < this.numY-1; j++)
                {
                    if (this.s[i, j] != 0)
                    {
                        u = (this.u[i, j] + this.u[i+1,j]) / 2;
                        v = (this.v[i, j] +
                             this.v[i, j + 1] ) / 2;
                            
                        x = i * this.height+ 0.5 * this.height - dt * u;
                        y = j * this.height + 0.5 * this.height - dt * v;
                        x = Math.Max(Math.Min(x, this.numX * height), height); //we have to be sure to not go out of the field
                        y = Math.Max(Math.Min(y, this.numY * height), height);
               

                        double x0 = this.height * 0.5;
                        double y0 = this.height * 0.5;

                        cell0x = (int)Math.Min(Math.Floor((x -x0) / height), this.numX - 1); //we calculate the position relative to the vectors but in cells to know the cell from wich we must get the velocity but taking care not to go out of the field
                        cell0y = (int)Math.Min(Math.Floor((y - y0) / height), this.numY - 1);
                           

                        int cell1x = Math.Min(cell0x + 1, this.numX - 1); //be sure the next cell is not out of the field either
                        int cell1y = Math.Min(cell0y + 1, this.numY - 1);

                        double b001 = ((x-x0) - cell0x * height) / height; //we calculate the position of each component inside the cell to make the weights
                        double b011 = ((y - y0 )- cell0y * height) / height;
                        double b000 = 1 - b001;
                        double b010 = 1 - b011;
                        auxSmokeField[i, j] = this.smokeField[cell0x, cell0y] * b000 * b010 +
                                     this.smokeField[cell1x, cell0y] * b001 * b010 +
                                     this.smokeField[cell0x, cell1y] * b000 * b011 +
                                     this.smokeField[cell1x, cell1y] * b001 * b011;
                        
                        /*Debug.Log("after: " +  smokeField[i,j,k]);*/

                    }
                }
            }
    
            this.smokeField = auxSmokeField;
           
        }        
    }
}