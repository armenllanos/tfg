using System;
using System.Collections.Generic;
using System.Threading;
using Unity.VisualScripting;
using UnityEngine;

namespace DefaultNamespace
{
    public class Fluid3D
    {
        private double height;
        private double density;
        private int numX;
        private int numY;
        private int numZ;
        private int numCells;
        double[,,] u;
        double[,,] v;
        double[,,] w;
        private float[] newU;
        private float[] newV;
        private double[,,] p;
        int[,,] s;
        public double[,,] smokeField;
        private float[] newM;
        private double overrelaxation = 1.9999;
        private string mode;
        private float speed;
        static double[,,] auxU;
        static double[,,] auxV;
        static double[,,] auxW;

        public Fluid3D(double density, int numX, int numY, int numZ, float height, float speed, string mode)
        {
            this.mode = mode;
            this.speed = speed;
            this.density = density;
            this.height = height;
            this.numX = numX + 2; //para añadir las celdas de los bordes
            this.numY = numY + 2;
            this.numZ = numZ + 2;
            this.numCells = this.numY * this.numX;
            this.u = new double[this.numX, this.numY, this.numZ];
            this.v = new double[this.numX, this.numY, this.numZ];
            this.w = new double[this.numX, this.numY, this.numZ];
            this.p = new double[this.numX, this.numY, this.numZ];
            this.s = new int[this.numX, this.numY, this.numZ];
            this.smokeField = new double[this.numX, this.numY, this.numZ];
            Fluid3D.auxU = new double[this.numX, this.numY, this.numZ];
            Fluid3D.auxV = new double[this.numX, this.numY, this.numZ];
            Fluid3D.auxW = new double[this.numX, this.numY, this.numZ];
            for (int i = 0; i < numX; i++) //we want to have limits, so we put 0 in the borders an 1 inside the field
            {
                for (int j = 0; j < numY; j++)
                {
                    for (int k = 0; k < numZ; k++)
                    {
                        /*if (k == 0 || k == this.numX-1)
                        {
                            this.s[i, j, k] = 0;
                        }
                        else if (j == 0 || j == this.numX-1)
                        {
                            this.s[i, j, k] = 0;
                        }
                        else if (i == 0 || i == this.numX-1)
                        {
                            this.s[i, j, k] = 0;
                        }
                        else
                        {
                            this.s[i, j, k] = 1;
                        }*/
                        this.s[i, j, k] = 0;
                        this.smokeField[i, j, k] = speed;
                        this.p[i, j, k] = 0;
                        this.u[i, j, k] = 0;
                        this.v[i, j, k] = 0;
                        this.w[i, j, k] = 0;
                    }
                }
            }

            for (int i = 1;
                 i < numX - 1;
                 i++) //we want to have limits, so we put 0 in the borders an 1 inside the field
            {
                for (int j = 1; j < numY - 1; j++)
                {
                    for (int k = 1; k < numZ - 1; k++)
                    {
                        /*if (k == 0 || k == this.numX-1)
                        {
                            this.s[i, j, k] = 0;
                        }
                        else if (j == 0 || j == this.numX-1)
                        {
                            this.s[i, j, k] = 0;
                        }
                        else if (i == 0 || i == this.numX-1)
                        {
                            this.s[i, j, k] = 0;
                        }
                        else
                        {
                            this.s[i, j, k] = 1;
                        }*/
                        this.s[i, j, k] = 1;
                    }
                }
            }


            //this.w = new float[this.numX,this.numY,this.numZ];
        }

        public void simulate(double dt)
        {
            if (mode.Equals("y0"))
            {
                for (int i = 0;
                     i < numX;
                     i++) //we want to have limits, so we put 0 in the borders an 1 inside the field
                {
                    for (int k = 0; k < numZ; k++)
                    {
                        v[i, 1, k] = speed;
                    }
                }
            }
            else if (mode.Equals("y1"))
            {
            }
            else if (mode.Equals("x0"))
            {
            }
            else if (mode.Equals("x1"))
            {
            }
            else if (mode.Equals("z0"))
            {
            }
            else if (mode.Equals("z1"))
            {
            }

            modifyVelocity(dt, -9.8);
            forceIncomprensibility(500, dt);
            advection(dt);
            smoke(dt);
        }

        public void modifyVelocity(double dt, double aceleration)
        {
            //first we need to take into account the base values for the different vectors, adding the gravity to the v vector
            for (int i = 1; i < this.numX - 1; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    for (int k = 1; k < this.numZ - 1; k++)
                    {
                        if (this.s[i, j, k] != 0.0 && this.s[i, j - 1, k] != 0.0)
                            this.v[i, j, k] += aceleration * dt;
                    }
                }
            }
        }

        class ThreadWorker
        {
            private int initiali;
            private int finali;
            private int initialj;
            private int finalj;
            private int initialk;
            private int finalk;
            private double dt;
            private Fluid3D fluid;

            public ThreadWorker(int initiali, int finali, int initialj, int finalj, int initialk, int finalk, double dt,
                Fluid3D fluid)
            {
                this.initiali = initiali;
                this.finali = finali;
                this.initialj = initialj;
                this.finalj = finalj;
                this.initialk = initialk;
                this.finalk = finalk;
                this.dt = dt;
                this.fluid = fluid;
            }

            /*public void ThreadGaussSiedel()
            {
                double outflow = 0;
                int actualS = 0;
                if (fluid.s[i, j, k] != 0)
                {
                    actualS = fluid.s[i + 1, j, k] + fluid.s[i - 1, j, k] + fluid.s[i, j + 1, k] +
                              fluid.s[i, j - 1, k] + fluid.s[i, j, k - 1] + fluid.s[i, j, k + 1];
                    if (actualS != 0)
                    {
                        outflow = fluid.u[i + 1, j, k] - fluid.u[i, j, k] + fluid.v[i, j + 1, k] -
                            fluid.v[i, j, k] + fluid.w[i, j, k + 1] - fluid.w[i, j, k];

                        fluid.u[i, j, k] = fluid.u[i, j, k] + outflow * fluid.s[i - 1, j, k] / actualS;
                        fluid.u[i + 1, j, k] =
                            fluid.u[i + 1, j, k] - outflow * fluid.s[i + 1, j, k] / actualS;
                        fluid.v[i, j, k] = fluid.v[i, j, k] + outflow * fluid.s[i, j - 1, k] / actualS;
                        fluid.v[i, j + 1, k] =
                            fluid.v[i, j + 1, k] - outflow * fluid.s[i, j + 1, k] / actualS;
                        fluid.w[i, j, k] = fluid.w[i, j, k] + outflow * fluid.s[i, j, k - 1] / actualS;
                        fluid.w[i, j, k + 1] =
                            fluid.w[i, j, k + 1] - outflow * fluid.s[i, j, k + 1] / actualS;

                        outflow *= fluid.overrelaxation;
                        fluid.p[i, j, k] = fluid.p[i, j, k] + outflow / actualS * fluid.density * fluid.height / dt;
                    }
                }
            }*/


            public void ThreadGaussSiedel()
            {
                double outflow = 0;
                int actualS = 0;

                for (int i = initiali; i < finali; i++)
                {
                    for (int j = initialj; j < finalj; j++)
                    {
                        for (int k = initialk; k < finalk; k++)
                        {
                            if (fluid.s[i, j, k] != 0)
                            {
                                actualS = fluid.s[i + 1, j, k] + fluid.s[i - 1, j, k] + fluid.s[i, j + 1, k] +
                                          fluid.s[i, j - 1, k] + fluid.s[i, j, k - 1] + fluid.s[i, j, k + 1];
                                if (actualS != 0)
                                {
                                    outflow = fluid.u[i + 1, j, k] - fluid.u[i, j, k] + fluid.v[i, j + 1, k] -
                                        fluid.v[i, j, k] + fluid.w[i, j, k + 1] - fluid.w[i, j, k];

                                    Fluid3D.auxU[i, j, k] =
                                        fluid.u[i, j, k] + outflow * fluid.s[i - 1, j, k] / actualS;
                                    /*converged = converged &&
                                                Math.Abs((auxU[i, j, k] - this.u[i, j, k]) / auxU[i, j, k]) < error;*/
                                    Fluid3D.auxU[i + 1, j, k] =
                                        fluid.u[i + 1, j, k] - outflow * fluid.s[i + 1, j, k] / actualS;
                                    /*converged = converged && Math.Abs((auxU[i+1, j, k] - this.u[i+1, j, k])/auxU[i+1, j, k]) < error;*/
                                    Fluid3D.auxV[i, j, k] =
                                        fluid.v[i, j, k] + outflow * fluid.s[i, j - 1, k] / actualS;
                                    /*converged = converged && Math.Abs((auxV[i, j, k] - this.v[i, j, k])/auxV[i, j, k]) < error;*/
                                    Fluid3D.auxV[i, j + 1, k] =
                                        fluid.v[i, j + 1, k] - outflow * fluid.s[i, j + 1, k] / actualS;
                                    /*converged = converged && Math.Abs((auxV[i, j+1, k] - this.v[i, j+1, k])/auxV[i, j+1, k]) < error;*/
                                    Fluid3D.auxW[i, j, k] =
                                        fluid.w[i, j, k] + outflow * fluid.s[i, j, k - 1] / actualS;
                                    /*converged = converged && Math.Abs((auxW[i, j, k] - this.w[i, j, k])/auxW[i, j, k]) < error;*/
                                    Fluid3D.auxW[i, j, k + 1] =
                                        fluid.w[i, j, k + 1] - outflow * fluid.s[i, j, k + 1] / actualS;
                                    /*converged = converged && Math.Abs((auxW[i, j, k+1] - this.w[i, j, k+1])/auxW[i, j, k+1]) < error;*/

                                    /*outflow *= overrelaxation;
                                    this.p[i, j, k] = p[i, j, k] + outflow / actualS * density * height / dt;*/
                                }
                            }
                        }
                    }
                }
            }
            public void ThreadGaussSiedelForU()
            {
                double outflow = 0;
                int actualS = 0;

                for (int i = 1; i < fluid.numX-1; i++)
                {
                    for (int j = 1; j < fluid.numY-1; j++)
                    {
                        for (int k = 1; k < fluid.numZ-1; k++)
                        {
                            if (fluid.s[i, j, k] != 0)
                            {
                                actualS = fluid.s[i + 1, j, k] + fluid.s[i - 1, j, k] + fluid.s[i, j + 1, k] +
                                          fluid.s[i, j - 1, k] + fluid.s[i, j, k - 1] + fluid.s[i, j, k + 1];
                                if (actualS != 0)
                                {
                                    outflow = fluid.u[i + 1, j, k] - fluid.u[i, j, k] + fluid.v[i, j + 1, k] -
                                        fluid.v[i, j, k] + fluid.w[i, j, k + 1] - fluid.w[i, j, k];

                                    Fluid3D.auxU[i, j, k] =
                                        fluid.u[i, j, k] + outflow * fluid.s[i - 1, j, k] / actualS;
                                    /*converged = converged &&
                                                Math.Abs((auxU[i, j, k] - this.u[i, j, k]) / auxU[i, j, k]) < error;*/
                                    Fluid3D.auxU[i + 1, j, k] =
                                        fluid.u[i + 1, j, k] - outflow * fluid.s[i + 1, j, k] / actualS;
                                    /*converged = converged && Math.Abs((auxW[i, j, k+1] - this.w[i, j, k+1])/auxW[i, j, k+1]) < error;*/

                                    /*outflow *= overrelaxation;
                                    this.p[i, j, k] = p[i, j, k] + outflow / actualS * density * height / dt;*/
                                }
                            }
                        }
                    }
                }
            }
            public void ThreadGaussSiedelForV()
            {
                double outflow = 0;
                int actualS = 0;

                for (int i = 1; i < fluid.numX-1; i++)
                {
                    for (int j = 1; j < fluid.numY-1; j++)
                    {
                        for (int k = 1; k < fluid.numZ-1; k++)
                        {
                            if (fluid.s[i, j, k] != 0)
                            {
                                actualS = fluid.s[i + 1, j, k] + fluid.s[i - 1, j, k] + fluid.s[i, j + 1, k] +
                                          fluid.s[i, j - 1, k] + fluid.s[i, j, k - 1] + fluid.s[i, j, k + 1];
                                if (actualS != 0)
                                {
                                    outflow = fluid.u[i + 1, j, k] - fluid.u[i, j, k] + fluid.v[i, j + 1, k] -
                                        fluid.v[i, j, k] + fluid.w[i, j, k + 1] - fluid.w[i, j, k];

                                    Fluid3D.auxV[i, j, k] =
                                        fluid.v[i, j, k] + outflow * fluid.s[i - 1, j, k] / actualS;
                                    /*converged = converged &&
                                                Math.Abs((auxU[i, j, k] - this.u[i, j, k]) / auxU[i, j, k]) < error;*/
                                    Fluid3D.auxV[i + 1, j, k] =
                                        fluid.v[i + 1, j, k] - outflow * fluid.s[i + 1, j, k] / actualS;
                                    /*converged = converged && Math.Abs((auxW[i, j, k+1] - this.w[i, j, k+1])/auxW[i, j, k+1]) < error;*/

                                    /*outflow *= overrelaxation;
                                    this.p[i, j, k] = p[i, j, k] + outflow / actualS * density * height / dt;*/
                                }
                            }
                        }
                    }
                }
            }
            public void ThreadGaussSiedelForW()
            {
                double outflow = 0;
                int actualS = 0;

                for (int i = 1; i < fluid.numX-1; i++)
                {
                    for (int j = 1; j < fluid.numY-1; j++)
                    {
                        for (int k = 1; k < fluid.numZ-1; k++)
                        {
                            if (fluid.s[i, j, k] != 0)
                            {
                                actualS = fluid.s[i + 1, j, k] + fluid.s[i - 1, j, k] + fluid.s[i, j + 1, k] +
                                          fluid.s[i, j - 1, k] + fluid.s[i, j, k - 1] + fluid.s[i, j, k + 1];
                                if (actualS != 0)
                                {
                                    outflow = fluid.u[i + 1, j, k] - fluid.u[i, j, k] + fluid.v[i, j + 1, k] -
                                        fluid.v[i, j, k] + fluid.w[i, j, k + 1] - fluid.w[i, j, k];

                                    Fluid3D.auxW[i, j, k] =
                                        fluid.w[i, j, k] + outflow * fluid.s[i - 1, j, k] / actualS;
                                    /*converged = converged &&
                                                Math.Abs((auxU[i, j, k] - this.u[i, j, k]) / auxU[i, j, k]) < error;*/
                                    Fluid3D.auxU[i + 1, j, k] =
                                        fluid.w[i + 1, j, k] - outflow * fluid.s[i + 1, j, k] / actualS;
                                    /*converged = converged && Math.Abs((auxW[i, j, k+1] - this.w[i, j, k+1])/auxW[i, j, k+1]) < error;*/

                                    /*outflow *= overrelaxation;
                                    this.p[i, j, k] = p[i, j, k] + outflow / actualS * density * height / dt;*/
                                }
                            }
                        }
                    }
                }
            }
        }


        /*public void forceIncomprensibilityMultyThread(int iterations, double dt)
        {
            for (int i = 1; i <= this.numX - 1; ++i)
            {
                for (int j = 1; j <= this.numY - 1; ++j)
                {
                    for (int k = 1; k <= this.numZ - 1; ++k)
                    {
                        if ((i + j + k) % 2 == 1)
                        {
                            ThreadWorker tw = new ThreadWorker(i, j, k, dt, this);
                            Thread t = new Thread(new ThreadStart(tw.ThreadGaussSiedel));
                            t.Start();
                            t.Join();
                        }
                    }
                }

                for (int j = 1; j <= this.numY - 1; ++j)
                {
                    for (int k = 1; k <= this.numZ - 1; ++k)
                    {
                        if ((i + j + k) % 2 == 0)
                        {
                            ThreadWorker tw = new ThreadWorker(i, j, k, dt, this);
                            Thread t = new Thread(new ThreadStart(tw.ThreadGaussSiedel));
                            t.Start();
                            t.Join();
                        }
                    }
                }
            }
        }*/

        public void forceIncomprensibility(int iterations, double dt)
        {
            double outflow = 0;
            int actualS = 0;
            double[,,] auxU = new double[this.numX, this.numY, this.numZ],
                auxV = new double[this.numX, this.numY, this.numZ],
                auxW = new double[this.numX, this.numY, this.numZ];
            bool converged = false;
            double error = 1;
            int l = 0;
          

            for (l = 0; (l < iterations ); l++)
            {
                /*List<Thread> threads = new List<Thread>();
                ThreadWorker tw1 = new ThreadWorker(1, this.numX / 2, 1, this.numY / 2,
                    1, this.numZ - 1, dt,
                    this);
                Thread t1 = new Thread(new ThreadStart(tw1.ThreadGaussSiedelForU));
                t1.Start();
                threads.Add(t1);
                ThreadWorker tw2 = new ThreadWorker(this.numX / 2, this.numX - 1, 1, this.numY / 2,
                    1, this.numZ - 1, dt,
                    this);
                Thread t2 = new Thread(new ThreadStart(tw2.ThreadGaussSiedelForV));
                t2.Start();
                threads.Add(t2);
                ThreadWorker tw3 = new ThreadWorker(1, this.numX / 2, this.numY / 2, this.numY - 1,
                    1, this.numZ - 1, dt, 
                    this);
                Thread t3 = new Thread(new ThreadStart(tw3.ThreadGaussSiedelForW));
                t3.Start();
                threads.Add(t3);
                */
                

                /*
                ThreadWorker tw1 = new ThreadWorker(1, this.numX / 2, 1, this.numY / 2,
                    1, this.numZ - 1, dt,
                    this);
                Thread t1 = new Thread(new ThreadStart(tw1.ThreadGaussSiedel));
                t1.Start();
                threads.Add(t1);
                ThreadWorker tw2 = new ThreadWorker(this.numX / 2, this.numX - 1, 1, this.numY / 2,
                    1, this.numZ - 1, dt,
                    this);
                Thread t2 = new Thread(new ThreadStart(tw2.ThreadGaussSiedel));
                t2.Start();
                threads.Add(t2);
                ThreadWorker tw3 = new ThreadWorker(1, this.numX / 2, this.numY / 2, this.numY - 1,
                    1, this.numZ - 1, dt, 
                    this);
                Thread t3 = new Thread(new ThreadStart(tw3.ThreadGaussSiedel));
                t3.Start();
                threads.Add(t3);
                ThreadWorker tw4 = new ThreadWorker(this.numX / 2, this.numX - 1, this.numY / 2, this.numY - 1,
                    1, this.numZ - 1, dt,
                    this);
                Thread t4 = new Thread(new ThreadStart(tw4.ThreadGaussSiedel));
                t4.Start();
                threads.Add(t4);

                */
                /*
                for (int i = 0; i < 3; i++)
                {
                    threads[i].Join();
                }
                */
               
                /*converged = true;*/
                for (int i = 1; i < this.numX - 1; i++)
                {
                    for (int j = 1; j < this.numY - 1; j++)
                    {
                        for (int k = 1; k < this.numZ - 1; k++)
                        {
                            if (this.s[i, j, k] != 0)
                            {
                                actualS = this.s[i + 1, j, k] + this.s[i - 1, j, k] + this.s[i, j + 1, k] +
                                          this.s[i, j - 1, k] + this.s[i, j, k - 1] + this.s[i, j, k + 1];
                                if (actualS != 0)
                                {
                                    outflow = this.u[i + 1, j, k] - this.u[i, j, k] + this.v[i, j + 1, k] -
                                        this.v[i, j, k] + this.w[i, j, k + 1] - this.w[i, j, k];

                                    auxU[i, j, k] = this.u[i, j, k] + outflow * this.s[i - 1, j, k] / actualS;
                                    /*converged = converged &&
                                                Math.Abs((auxU[i, j, k] - this.u[i, j, k]) / auxU[i, j, k]) < error;*/
                                    auxU[i + 1, j, k] =
                                        this.u[i + 1, j, k] - outflow * this.s[i + 1, j, k] / actualS;
                                    /*converged = converged && Math.Abs((auxU[i+1, j, k] - this.u[i+1, j, k])/auxU[i+1, j, k]) < error;*/
                                    auxV[i, j, k] = this.v[i, j, k] + outflow * this.s[i, j - 1, k] / actualS;
                                    /*converged = converged && Math.Abs((auxV[i, j, k] - this.v[i, j, k])/auxV[i, j, k]) < error;*/
                                    auxV[i, j + 1, k] =
                                        this.v[i, j + 1, k] - outflow * this.s[i, j + 1, k] / actualS;
                                    /*converged = converged && Math.Abs((auxV[i, j+1, k] - this.v[i, j+1, k])/auxV[i, j+1, k]) < error;*/
                                    auxW[i, j, k] = this.w[i, j, k] + outflow * this.s[i, j, k - 1] / actualS;
                                    /*converged = converged && Math.Abs((auxW[i, j, k] - this.w[i, j, k])/auxW[i, j, k]) < error;*/
                                    auxW[i, j, k + 1] =
                                        this.w[i, j, k + 1] - outflow * this.s[i, j, k + 1] / actualS;
                                    /*converged = converged && Math.Abs((auxW[i, j, k+1] - this.w[i, j, k+1])/auxW[i, j, k+1]) < error;*/

                                    /*outflow *= overrelaxation;
                                    this.p[i, j, k] = p[i, j, k] + outflow / actualS * density * height / dt;*/
                                }
                            }
                        }
                    }
                }
                u = auxU;
                v = auxV;
                w = auxW;
              
            }
            Debug.Log("n iteraciones = " + l);

            for (int k = 0; k < this.numZ - 1; k++)
            {
                if (k == 0)
                {
                    for (int i = 1; i < this.numX - 1; i++)
                    {
                        for (var j = 1; j < this.numY - 1; j++)
                        {
                            this.w[i, j, k] = this.w[i, j, k + 1];
                        }
                    }
                }
                else if (k == this.numZ - 1)
                {
                    for (int i = 1; i < this.numX - 1; i++)
                    {
                        for (var j = 1; j < this.numY - 1; j++)
                        {
                            this.w[i, j, k] = this.w[i, j, k - 1];
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < this.numX; i++)
                    {
                        this.u[i, 0, 0] = this.u[i, 1, k];
                        this.u[i, this.numY - 1, k] = this.u[i, this.numY - 2, k];
                    }

                    for (var j = 0; j < this.numY; j++)
                    {
                        this.v[0, j, k] = this.v[1, j, k];
                        this.v[this.numX - 1, j, k] = this.v[this.numX - 2, j, k];
                    }
                }
            }
        }

        public void advection(double dt)
        {
            double x, y, z, x0, y0, z0;
            double u, v, w;
            int u0, v0, w0;
            int i0, j0, k0, i1, j1, k1;
            double velocity;
            double[,,] auxu = new double[this.numX, this.numY, this.numZ];
            double[,,] auxv = new double[this.numX, this.numY, this.numZ];
            double[,,] auxw = new double[this.numX, this.numY, this.numZ];
            double b000, b001, b010, b011, b100, b101;
            for (int i = 1; i < this.numX - 1; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    for (int k = 1; k < this.numZ - 1; k++)
                    {
                        //u
                        if (this.s[i, j, k] != 0.0 && this.s[i - 1, j, k] != 0.0)
                        {
                            u = this.u[i, j, k];
                            v = (this.v[i, j, k] + this.v[i - 1, j, k] +
                                 this.v[i, j + 1, k] + this.v[i - 1, j + 1, k]) / 4;
                            w = (this.w[i, j, k] + this.w[i - 1, j, k] +
                                 this.w[i, j, k + 1] + this.w[i - 1, j, k + 1]) / 4;

                            x = i * this.height - dt * u;
                            z = k * this.height + 0.5 * this.height - dt * w;
                            y = j * this.height + 0.5 * this.height - dt * v;
                            x = Math.Max(Math.Min(x, this.numX * height),
                                height); //we have to be sure to not go out of the field
                            y = Math.Max(Math.Min(y, this.numY * height), height);
                            z = Math.Max(Math.Min(z, this.numZ * height), height);

                            x0 = 0;
                            z0 = this.height * 0.5;
                            y0 = this.height * 0.5;

                            i0 = (int)Math.Min(Math.Floor(x / height),
                                this.numX -
                                1); //we calculate the position relative to the vectors but in cells to know the cell from wich we must get the velocity but taking care not to go out of the field
                            j0 = (int)Math.Min(Math.Floor((y - y0) / height), this.numY - 1);
                            k0 = (int)Math.Min(Math.Floor((z - z0) / height), this.numZ - 1);

                            i1 = Math.Min(i0 + 1,
                                this.numX - 1); //be sure the next cell is not out of the field either
                            j1 = Math.Min(j0 + 1, this.numY - 1);
                            k1 = Math.Min(k0 + 1, this.numZ - 1);

                            b001 = (x - i0 * height) /
                                   height; //we calculate the position of each component inside the cell to make the weights
                            b011 = (y - y0 - j0 * height) / height;
                            b101 = (z - z0 - k0 * height) / height;

                            b000 = 1 - b001;
                            b010 = 1 - b011;
                            b100 = 1 - b101;

                            auxu[i, j, k] = this.u[i0, j0, k0] * b000 * b100 * b010 +
                                            this.u[i1, j0, k0] * b001 * b100 * b010 +
                                            this.u[i1, j0, k1] * b001 * b101 * b010 +
                                            this.u[i1, j0, k1] * b000 * b101 * b010
                                            + this.u[i0, j1, k0] * b000 * b100 * b011 +
                                            this.u[i1, j1, k0] * b001 * b100 * b011
                                            + this.u[i0, j1, k1] * b000 * b101 * b011 +
                                            this.u[i1, j1, k1] * b001 * b101 * b011;
                        }


                        //v
                        if (this.s[i, j, k] != 0.0 && this.s[i, j - 1, k] != 0.0)
                        {
                            v = this.v[i, j, k];
                            u = (this.u[i, j - 1, k] + this.u[i, j, k] +
                                 this.u[i + 1, j - 1, k] + this.u[i + 1, j, k]) / 4;
                            w = (this.w[i, j - 1, k] + this.w[i, j, k] +
                                 this.w[i, j - 1, k + 1] + this.w[i, j, k + 1]) / 4;

                            x = i * this.height + 0.5 * this.height - dt * u;
                            z = k * this.height + 0.5 * this.height - dt * w;
                            y = j * this.height - dt * v;
                            x = Math.Max(Math.Min(x, this.numX * height),
                                height); //we have to be sure to not go out of the field
                            y = Math.Max(Math.Min(y, this.numY * height), height);
                            z = Math.Max(Math.Min(z, this.numZ * height), height);

                            x0 = this.height * 0.5;
                            z0 = this.height * 0.5;
                            y0 = 0;

                            i0 = (int)Math.Min(Math.Floor((x - x0) / height),
                                this.numX -
                                1); //we calculate the position relative to the vectors but in cells to know the cell from wich we must get the velocity but taking care not to go out of the field
                            j0 = (int)Math.Min(Math.Floor(y / height), this.numY - 1);
                            k0 = (int)Math.Min(Math.Floor((z - z0) / height), this.numZ - 1);

                            i1 = Math.Min(i0 + 1,
                                this.numX - 1); //be sure the next cell is not out of the field either
                            j1 = Math.Min(j0 + 1, this.numY - 1);
                            k1 = Math.Min(k0 + 1, this.numZ - 1);

                            b001 = (x - x0 - i0 * height) /
                                   height; //we calculate the position of each component inside the cell to make the weights
                            b011 = (y - j0 * height) / height;
                            b101 = (z - z0 - k0 * height) / height;

                            b000 = 1 - b001;
                            b010 = 1 - b011;
                            b100 = 1 - b101;

                            auxv[i, j, k] = this.v[i0, j0, k0] * b000 * b100 * b010 +
                                            this.v[i1, j0, k0] * b001 * b100 * b010 +
                                            this.v[i1, j0, k1] * b001 * b101 * b010 +
                                            this.v[i1, j0, k1] * b000 * b101 * b010
                                            + this.v[i0, j1, k0] * b000 * b100 * b011 +
                                            this.v[i1, j1, k0] * b001 * b100 * b011
                                            + this.v[i0, j1, k1] * b000 * b101 * b011 +
                                            this.v[i1, j1, k1] * b001 * b101 * b011;
                        }

                        //w
                        if (this.s[i, j, k] != 0.0 && this.s[i, j, k - 1] != 0.0)
                        {
                            v = (this.v[i, j, k] + this.v[i, j, k - 1] +
                                 this.v[i, j + 1, k] + this.v[i, j + 1, k - 1]) / 4;

                            u = (this.u[i, j, k - 1] + this.u[i, j, k] +
                                 this.u[i + 1, j, k - 1] + this.u[i + 1, j, k]) / 4;
                            w = this.w[i, j, k];

                            x = i * this.height + 0.5 * this.height - dt * u;
                            z = k * this.height - dt * w;
                            y = j * this.height + 0.5 * this.height - dt * v;
                            x = Math.Max(Math.Min(x, this.numX * height),
                                height); //we have to be sure to not go out of the field
                            y = Math.Max(Math.Min(y, this.numY * height), height);
                            z = Math.Max(Math.Min(z, this.numZ * height), height);

                            x0 = this.height * 0.5;
                            z0 = 0;
                            y0 = this.height * 0.5;

                            i0 = (int)Math.Min(Math.Floor((x - x0) / height),
                                this.numX -
                                1); //we calculate the position relative to the vectors but in cells to know the cell from wich we must get the velocity but taking care not to go out of the field
                            j0 = (int)Math.Min(Math.Floor((y - y0) / height), this.numY - 1);
                            k0 = (int)Math.Min(Math.Floor((z - z0) / height), this.numZ - 1);

                            i1 = Math.Min(i0 + 1,
                                this.numX - 1); //be sure the next cell is not out of the field either
                            j1 = Math.Min(j0 + 1, this.numY - 1);
                            k1 = Math.Min(k0 + 1, this.numZ - 1);

                            b001 = (x - x0 - i0 * height) /
                                   height; //we calculate the position of each component inside the cell to make the weights
                            b011 = (y - y0 - j0 * height) / height;
                            b101 = (z - k0 * height) / height;

                            b000 = 1 - b001;
                            b010 = 1 - b011;
                            b100 = 1 - b101;

                            auxw[i, j, k] = this.w[i0, j0, k0] * b000 * b100 * b010 +
                                            this.w[i1, j0, k0] * b001 * b100 * b010 +
                                            this.w[i1, j0, k1] * b001 * b101 * b010 +
                                            this.w[i1, j0, k1] * b000 * b101 * b010
                                            + this.w[i0, j1, k0] * b000 * b100 * b011 +
                                            this.w[i1, j1, k0] * b001 * b100 * b011
                                            + this.w[i0, j1, k1] * b000 * b101 * b011 +
                                            this.w[i1, j1, k1] * b001 * b101 * b011;
                        }


                        /*u0 = Math.Min((int)Math.Floor(x / this.height),1);
                        if (y % this.height > this.height/2)
                        {
                            if (z % this.height > this.height / 2)
                            {
                                x0 = x % this.height;
                                z0 = z % this.height - this.height * 0.5;
                                y0 = y % this.height - this.height * 0.5;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                j0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxu[i,j,k] = this.v[u0,j0,k0] * b000 * b100 * b010 + this.v[u0+1,j0,k0] * b001 * b100 * b010 +
                                              this.v[u0+1,j0,k0+1] * b001 * b101 * b010  + this.v[u0,j0,k0+1] * b000 * b101 * b010 
                                              + this.v[u0,j0+1,k0] * b000 * b100 * b011 + this.v[u0+1,j0+1,k0] * b001 * b100 * b011
                                              + this.v[u0,j0+1,k0+1] * b000 * b101 * b011 + this.v[u0+1,j0+1,k0+1] * b001 * b101 * b011;
                            }
                            else
                            {
                                x0 = x % this.height;
                                z0 = z % this.height;
                                y0 = y % this.height - this.height * 0.5;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                j0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                if (k0 == 0)
                                {
                                    k0++;
                                }
                                auxu[i,j,k] = this.v[u0,j0,k0-1] * b000 * b100 * b010 + this.v[u0+1,j0,k0-1] * b001 * b100 * b010 +
                                              this.v[u0+1,j0,k0] * b001 * b101 * b010  + this.v[u0,j0,k0] * b000 * b101 * b010 
                                              + this.v[u0,j0+1,k0-1] * b000 * b100 * b011 + this.v[u0+1,j0+1,k0-1] * b001 * b100 * b011
                                              + this.v[u0,j0+1,k0] * b000 * b101 * b011 + this.v[u0+1,j0+1,k0] * b001 * b101 * b011;
                                
                            }
                        }
                        else
                        {
                            if (z % this.height > this.height / 2)
                            {
                                x0 = x % this.height;
                                z0 = z % this.height - this.height * 0.5;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                j0 =  (int)Math.Floor(x/ this.height);
                                if (j0 == 0)
                                {
                                    j0++;
                                }
                                k0 = (int)Math.Floor(z/ this.height);
                                auxu[i,j,k] = this.v[u0,j0-1,k0] * b000 * b100 * b010 + this.v[u0+1,j0-1,k0] * b001 * b100 * b010 +
                                              this.v[u0+1,j0,k0+1] * b001 * b101 * b010  + this.v[u0,j0-1,k0+1] * b000 * b101 * b010 
                                              + this.v[u0,j0,k0] * b000 * b100 * b011 + this.v[u0+1,j0,k0] * b001 * b100 * b011
                                              + this.v[u0,j0,k0+1] * b000 * b101 * b011 + this.v[u0+1,j0,k0+1] * b001 * b101 * b011;
                            }
                            else
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                               
                                j0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z / this.height);
                                if (j0 == 0)
                                {
                                    j0++;
                                }
                                if (k0 == 0)
                                {
                                    k0++;
                                }
                                Debug.Log("x0: "+u0+" y0: "+j0+" z0: "+k0);
                                auxu[i,j,k] = this.v[u0,j0-1,k0-1] * b000 * b100 * b010 + this.v[u0+1,j0-1,k0-1] * b001 * b100 * b010 +
                                              this.v[u0+1,j0,k0] * b001 * b101 * b010  + this.v[u0,j0-1,k0] * b000 * b101 * b010 
                                              + this.v[u0,j0,k0-1] * b000 * b100 * b011 + this.v[u0+1,j0,k0-1] * b001 * b100 * b011
                                              + this.v[u0,j0,k0] * b000 * b101 * b011 + this.v[u0+1,j0,k0] * b001 * b101 * b011;
                                
                            }
                        }
                        //v
                        v = this.v[i, j,k];
                        u = (this.u[i, j - 1,k] + this.u[i, j,k] +
                             this.u[i + 1, j - 1,k] + this.u[i +1 , j,k]) / 4;
                        w = (this.w[i, j - 1,k] + this.w[i, j,k] +
                             this.w[i , j - 1,k+1] + this.w[i, j,k+1]) / 4;
                        x = i * this.height + 0.5 * this.height - dt * u;
                        z = k * this.height + 0.5 * this.height - dt * w;
                        y = j * this.height - dt * v;
                        v0 = (int)Math.Floor(y/ this.height);
                        if (x % this.height > this.height/2)
                        {
                            if (z % this.height > this.height / 2)
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height - this.height * 0.5;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxv[i,j,k] = this.v[i0,v0,k0] * b000 * b100 * b010 + this.v[i0+1,v0,k0] * b001 * b100 * b010 +
                                              this.v[i0+1,v0,k0+1] * b001 * b101 * b010  + this.v[i0,v0,k0+1] * b000 * b101 * b010 
                                              + this.v[i0,v0+1,k0] * b000 * b100 * b011 + this.v[i0+1,v0+1,k0] * b001 * b100 * b011
                                              + this.v[i0,v0+1,k0+1] * b000 * b101 * b011 + this.v[i0+1,v0+1,k0+1] * b001 * b101 * b011;
                            }
                            else
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxv[i,j,k] = this.v[i0,v0,k0-1] * b000 * b100 * b010 + this.v[i0+1,v0,k0-1] * b001 * b100 * b010 +
                                              this.v[i0+1,v0,k0] * b001 * b101 * b010  + this.v[i0,v0,k0] * b000 * b101 * b010 
                                              + this.v[i0,v0+1,k0-1] * b000 * b100 * b011 + this.v[i0+1,v0+1,k0-1] * b001 * b100 * b011
                                              + this.v[i0, v0 + 1, k0] * b000 * b101 * b011 +
                                              this.v[i0 + 1, v0 + 1, k0] * b001 * b101 * b011;
                                
                            }
                        }
                        else
                        {
                            if (z % this.height > this.height / 2)
                            {
                                x0 = x % this.height;
                                z0 = z % this.height - this.height * 0.5;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxv[i,j,k] = this.v[i0-1,v0,k0] * b000 * b100 * b010 + this.v[i0,v0,k0] * b001 * b100 * b010 +
                                              this.v[i0,v0,k0+1] * b001 * b101 * b010  + this.v[i0-1,v0,k0+1] * b000 * b101 * b010 
                                              + this.v[i0-1,v0+1,k0] * b000 * b100 * b011 + this.v[i0,v0+1,k0] * b001 * b100 * b011
                                              + this.v[i0-1,v0+1,k0+1] * b000 * b101 * b011 + this.v[i0,v0+1,k0+1] * b001 * b101 * b011;
                            }
                            else
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxv[i,j,k] = this.v[i0-1,v0,k0-1] * b000 * b100 * b010 + this.v[i0,v0,k0-1] * b001 * b100 * b010 +
                                              this.v[i0,v0,k0] * b001 * b101 * b010  + this.v[i0-1,v0,k0] * b000 * b101 * b010 
                                              + this.v[i0-1,v0+1,k0-1] * b000 * b100 * b011 + this.v[i0,v0+1,k0-1] * b001 * b100 * b011
                                              + this.v[i0-1, v0 + 1, k0] * b000 * b101 * b011 +
                                              this.v[i0, v0 + 1, k0] * b001 * b101 * b011;
                                
                            }
                        }
                        //w
                        v = (this.v[i, j ,k] + this.v[i, j,k-1] +
                             this.v[i, j + 1,k] + this.v[i, j+1,k-1]) / 4;;
                        u = (this.u[i, j,k-1] + this.u[i, j,k] +
                             this.u[i + 1, j ,k -1 ] + this.u[i +1 , j,k]) / 4;
                        w = this.w[i,j,k];
                        x = i * this.height + 0.5 * this.height - dt * u;
                        z = k * this.height  - dt * w;
                        y = j * this.height + 0.5 * this.height - dt * v;
                        w0 = (int)Math.Floor(z/ this.height);
                        if (x % this.height > this.height/2)
                        {
                            if (y % this.height > this.height / 2)
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height ;
                                y0 = y % this.height - this.height * 0.5;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                j0 = (int)Math.Floor(y/ this.height);
                                auxw[i,j,k] = this.w[i0,j0,w0] * b000 * b100 * b010 + this.w[i0+1,j0,w0] * b001 * b100 * b010 +
                                              this.w[i0+1,j0,w0+1] * b001 * b101 * b010  + this.w[i0,j0,w0+1] * b000 * b101 * b010 
                                              + this.w[i0,j0+1,w0] * b000 * b100 * b011 + this.w[i0+1,j0+1,w0] * b001 * b100 * b011
                                              + this.w[i0,j0+1,w0+1] * b000 * b101 * b011 + this.w[i0+1,j0+1,w0+1] * b001 * b101 * b011;
                            }
                            else
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                j0 = (int)Math.Floor(y/ this.height);
                                auxw[i,j,k] = this.w[i0,j0-1,w0] * b000 * b100 * b010 + this.w[i0+1,j0-1,w0] * b001 * b100 * b010 +
                                              this.w[i0+1,j0-1,w0+1] * b001 * b101 * b010  + this.w[i0,j0-1,w0+1] * b000 * b101 * b010 
                                              + this.w[i0,j0,w0] * b000 * b100 * b011 + this.w[i0+1,j0,w0] * b001 * b100 * b011
                                              + this.w[i0,j0,w0+1] * b000 * b101 * b011 + this.w[i0+1,j0,w0+1] * b001 * b101 * b011;
                            }
                        }
                        else
                        {
                            if (y % this.height > this.height / 2)
                            {
                                x0 = x % this.height;
                                z0 = z % this.height - this.height * 0.5;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                j0 = (int)Math.Floor(y/ this.height);
                                auxw[i,j,k] = this.w[i0-1,j0,w0] * b000 * b100 * b010 + this.w[i0,j0,w0] * b001 * b100 * b010 +
                                              this.w[i0,j0,w0+1] * b001 * b101 * b010  + this.w[i0-1,j0,w0+1] * b000 * b101 * b010 
                                              + this.w[i0-1,j0+1,w0] * b000 * b100 * b011 + this.w[i0,j0+1,w0] * b001 * b100 * b011
                                              + this.w[i0-1,j0+1,w0+1] * b000 * b101 * b011 + this.w[i0,j0+1,w0+1] * b001 * b101 * b011;
                            }
                            else
                            {
                                x0 = x % this.height;
                                z0 = z % this.height - this.height * 0.5;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                j0 = (int)Math.Floor(y/ this.height);
                                auxw[i,j,k] = this.w[i0-1,j0-1,w0] * b000 * b100 * b010 + this.w[i0,j0-1,w0] * b001 * b100 * b010 +
                                              this.w[i0,j0-1,w0+1] * b001 * b101 * b010  + this.w[i0-1,j0-1,w0+1] * b000 * b101 * b010 
                                              + this.w[i0-1,j0,w0] * b000 * b100 * b011 + this.w[i0,j0,w0] * b001 * b100 * b011
                                              + this.w[i0-1,j0,w0+1] * b000 * b101 * b011 + this.w[i0,j0,w0+1] * b001 * b101 * b011;
                                
                            }
                        }
                       
                        //v
                        v = this.v[i, j,k];
                        u = (this.u[i, j - 1,k] + this.u[i, j,k] +
                             this.u[i + 1, j - 1,k] + this.u[i +1 , j,k]) / 4;
                        w = (this.w[i, j - 1,k] + this.w[i, j,k] +
                             this.w[i , j - 1,k+1] + this.w[i, j,k+1]) / 4;
                        x = i * this.height + 0.5 * this.height - dt * u;
                        z = k * this.height + 0.5 * this.height - dt * w;
                        y = j * this.height - dt * v;
                        v0 = (int)Math.Floor(y/ this.height);
                        if (x % this.height > this.height/2)
                        {
                            if (z % this.height > this.height / 2)
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height - this.height * 0.5;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxv[i,j,k] = this.v[i0,v0,k0] * b000 * b100 * b010 + this.v[i0+1,v0,k0] * b001 * b100 * b010 +
                                              this.v[i0+1,v0,k0+1] * b001 * b101 * b010  + this.v[i0,v0,k0+1] * b000 * b101 * b010 
                                              + this.v[i0,v0+1,k0] * b000 * b100 * b011 + this.v[i0+1,v0+1,k0] * b001 * b100 * b011
                                              + this.v[i0,v0+1,k0+1] * b000 * b101 * b011 + this.v[i0+1,v0+1,k0+1] * b001 * b101 * b011;
                            }
                            else
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxv[i,j,k] = this.v[i0,v0,k0-1] * b000 * b100 * b010 + this.v[i0+1,v0,k0-1] * b001 * b100 * b010 +
                                              this.v[i0+1,v0,k0] * b001 * b101 * b010  + this.v[i0,v0,k0] * b000 * b101 * b010 
                                              + this.v[i0,v0+1,k0-1] * b000 * b100 * b011 + this.v[i0+1,v0+1,k0-1] * b001 * b100 * b011
                                              + this.v[i0, v0 + 1, k0] * b000 * b101 * b011 +
                                              this.v[i0 + 1, v0 + 1, k0] * b001 * b101 * b011;
                                
                            }
                        }
                        else
                        {
                            if (z % this.height > this.height / 2)
                            {
                                x0 = x % this.height;
                                z0 = z % this.height - this.height * 0.5;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxv[i,j,k] = this.v[i0-1,v0,k0] * b000 * b100 * b010 + this.v[i0,v0,k0] * b001 * b100 * b010 +
                                              this.v[i0,v0,k0+1] * b001 * b101 * b010  + this.v[i0-1,v0,k0+1] * b000 * b101 * b010 
                                              + this.v[i0-1,v0+1,k0] * b000 * b100 * b011 + this.v[i0,v0+1,k0] * b001 * b100 * b011
                                              + this.v[i0-1,v0+1,k0+1] * b000 * b101 * b011 + this.v[i0,v0+1,k0+1] * b001 * b101 * b011;
                            }
                            else
                            {
                                x0 = x % this.height - this.height * 0.5;
                                z0 = z % this.height;
                                y0 = y % this.height;
                                b000 = 1 - x0 / this.height;
                                b001 = x0 / this.height;
                                b010 = 1 - y0 / this.height;
                                b011 = y0 / this.height;
                                b100 = 1 - z0 / this.height;
                                b101 = z0 / this.height;
                                i0 =  (int)Math.Floor(x/ this.height);
                                k0 = (int)Math.Floor(z/ this.height);
                                auxv[i,j,k] = this.v[i0-1,v0,k0-1] * b000 * b100 * b010 + this.v[i0,v0,k0-1] * b001 * b100 * b010 +
                                              this.v[i0,v0,k0] * b001 * b101 * b010  + this.v[i0-1,v0,k0] * b000 * b101 * b010 
                                              + this.v[i0-1,v0+1,k0-1] * b000 * b100 * b011 + this.v[i0,v0+1,k0-1] * b001 * b100 * b011
                                              + this.v[i0-1, v0 + 1, k0] * b000 * b101 * b011 +
                                              this.v[i0, v0 + 1, k0] * b001 * b101 * b011;
                                
                            }
                        }*/
                    }
                }
            }

            this.v = auxv;
            this.u = auxu;
            this.w = auxw;
        }

        public void smoke(double dt)
        {
            double u, v, w;
            double x, y, z;
            int cell0x;
            int cell0y;
            int cell0z;
            double[,,] auxSmokeField = new double[this.numX, this.numY, this.numZ];
            for (int i = 1; i < this.numX - 1; i++)
            {
                for (int j = 1; j < this.numY - 1; j++)
                {
                    for (int k = 1; k < this.numZ - 1; k++)
                    {
                        if (this.s[i, j, k] != 0)
                        {
                            u = (this.u[i, j, k] + this.u[i + 1, j, k]) * 0.5;
                            v = (this.v[i, j, k] + this.v[i, j + 1, k]) * 0.5;
                            w = (this.w[i, j, k] + this.w[i, j, k + 1]) * 0.5;
                            x = i * this.height + this.height * 0.5 - dt * u;
                            y = j * this.height + this.height * 0.5 - dt * v;
                            z = k * this.height + this.height * 0.5 - dt * w;
                            cell0x = (int)Math.Floor((x - this.height * 0.5) / this.height);
                            cell0y = (int)Math.Floor((y - this.height * 0.5) / this.height);
                            cell0z = (int)Math.Floor((z - this.height * 0.5) / this.height);

                            cell0x = (int)Math.Min(Math.Max(0, cell0x), this.numX - 1);
                            cell0y = (int)Math.Min(Math.Max(0, cell0y), this.numY - 1);
                            cell0z = (int)Math.Min(Math.Max(0, cell0z), this.numZ - 1);

                            double b001 = ((x - this.height * 0.5) - cell0x * this.height) / this.height; //peso 
                            double b000 = 1 - b001;

                            double b011 = ((y - this.height * 0.5) - cell0y * this.height) / this.height;
                            double b010 = 1 - b011;
                            double b101 = ((z - this.height * 0.5) - cell0z * this.height) / this.height;
                            double b100 = 1 - b101;
                            int cell1y = Math.Min(cell0y + 1, this.numY - 1);
                            int cell1x = Math.Min(cell0x + 1, this.numX - 1);
                            int cell1z = Math.Min(cell0x + 1, this.numZ - 1);
                            /*Debug.Log("x: " + cell0x + "," + cell1x + "y: " + cell0y + "," + cell1y + "z: " +
                                      cell0z +
                                      "," + cell1z);*/
                            auxSmokeField[i, j, k] = b000 * b100 * b010 * this.smokeField[cell0x, cell0y, cell0z] +
                                                     b001 * b010 * b100 * this.smokeField[cell1x, cell0y, cell0z]
                                                     + b001 * b010 * b101 *
                                                     this.smokeField[cell1x, cell0y, cell1z] +
                                                     b000 * b010 * b101 * this.smokeField[cell0x, cell0y, cell1z]
                                                     + b000 * b011 * b100 *
                                                     this.smokeField[cell0x, cell1y, cell0z] +
                                                     b001 * b011 * b100 * this.smokeField[cell1x, cell1y, cell0z]
                                                     + b001 * b011 * b101 *
                                                     this.smokeField[cell0x, cell1y, cell1z] +
                                                     b001 * b011 * b101 * this.smokeField[cell1x, cell1y, cell1z];
                        }
                    }
                }
            }

            this.smokeField = auxSmokeField;
        }
    }
}