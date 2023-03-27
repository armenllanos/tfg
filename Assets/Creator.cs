using System;
using System.Collections;
using System.Collections.Generic;
using DefaultNamespace;
using UnityEngine;

public class Creator : MonoBehaviour
{
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
    float opacity = 0;

    void Start()
    {
        fluid = new Fluid3D(1000, numCubesX, numCubesY, numCubesZ, 1, speed, "y0");
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
        fluid.simulate(Time.deltaTime);
        for (int x = 1; x < numCubesX - 1; x++)
        {
            for (int y = 1; y < numCubesY - 1; y++)
            {
                for (int z = 1; z < numCubesZ - 1; z++)
                {
                    myColor = cubes[x, y, z].GetComponent<Renderer>().material.color;
                    if ((float)fluid.smokeField[x, y, z] <= speed)
                    {
                        cubes[x, y, z].GetComponent<Renderer>().material.color =
                            new Color(myColor.r, myColor.g, myColor.b, 1 - (float)fluid.smokeField[x, y, z] / speed);
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
}