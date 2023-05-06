using System;
using System.Collections;
using System.Collections.Generic;
using DefaultNamespace;using OpenCover.Framework.Model;
using UnityEngine;

public class Creator : MonoBehaviour
{

    public ComputeShader computeShader;

    public RenderTexture renderTexture;
    // Start is called before the first frame update
    public int numCubesX = 10;
    public int numCubesZ = 10;
    public int numCubesY = 10;
    public float spacing = 1f;
    public GameObject cubePrefab;
    public float speed;
    private GameObject[,,] cubes;
    private Fluid3D fluid;
    private ArrayList gameObjects = new ArrayList();
    private int j = 0;
    public bool GPU = false;
    float opacity = 0;

    private double[,,] smoke;

    void Start()
    {
        fluid = new Fluid3D(1000, numCubesX, numCubesY, numCubesZ, 1, speed, "y0",computeShader);
        cubes = new GameObject[numCubesX, numCubesY, numCubesZ];
        for (int x = 0; x < numCubesX; x++)
        {
            for (int y = 0; y < numCubesY; y++)
            {
                for (int z = 0; z < numCubesZ; z++)
                {
                    
                    Vector3 position = new Vector3(x * spacing + (float)0.5, y * spacing + (float)0.5,
                        z * spacing + (float)0.5);
                    cubes[x, y, z] = Instantiate(cubePrefab, position, Quaternion.identity, transform);
                    gameObjects.Add(cubes[x, y, z]);
                }
            }
        }
    }

    // Update is called once per frame
    void Update()
    {
        Color myColor;
        if (GPU)
        {
            fluid.simulateGPU(Time.deltaTime);
        }
        else
        {
            fluid.simulate(Time.deltaTime);
        }
        
        for (int x = 0; x < numCubesX; x++)
        {
            for (int y = 0; y < numCubesY; y++)
            {
                for (int z = 0; z < numCubesZ; z++)
                {
                    myColor = cubes[x, y, z].GetComponent<Renderer>().material.color;
                    Color finalColor = getColor(fluid.smokeField[x+1,y+1,z+1]);
                    /*Debug.Log("opacity: " + fluid.smokeField[x, y, z]);*/
              
                    
                    if ((float)fluid.smokeField[x+1, y+1, z+1] <= speed)
                    {
                        cubes[x, y, z].GetComponent<Renderer>().material.color =
                            finalColor;
                    }
                    else
                    {
                        cubes[x, y, z].GetComponent<Renderer>().material.color =
                            new Color(myColor.r, myColor.g, myColor.b, 0);
                    }

                    /*Debug.Log("opacity: " + fluid.smokeField[x, y, z] / speed + "smoke value" +
                             fluid.smokeField[x, y, z] / speed);*/
                }
            }
        }
    }

    private Color getColor(double val)
    {
        val = Math.Min(Math.Max(val, 0), 1- 0.0001);
        int d = 1;
        //val = d == 0.0 ? 0.5 : (val - 0) / d;
        double m = 0.25;
        int num = (int)Math.Floor(val / m);
        double s = (val - num * m) / m;
        /*double r=0, g=0, b=0;

        switch (num) {
            case 0 : r = 0.0; g = s; b = 1.0; break;
            case 1 : r = 0.0; g = 1.0; b = 1.0-s; break;
            case 2 : r = s; g = 1.0; b = 0.0; break;
            case 3 : r = 1.0; g = 1.0 - s; b = 0.0; break;
        }

        return new Color((float)(255 * r), (float)(255 * g), (float)(255 * b),(float)0.5);*/
        return new Color((float)0.5, (float)0.5, (float)0.5,(float)((1-val)/1.8));
    }
}
