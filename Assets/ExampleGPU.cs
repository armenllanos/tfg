using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ExampleGPU : MonoBehaviour
{
    // Start is called before the first frame update
    public ComputeShader computeShader;
    
    void Start()
    {
        int size = 20;
        int sizey = 10;
        int sizez = 10;
        int[] matriz = new int[size* sizey* sizez];
        int[] matriz1 = new int[size* sizey* sizez];
        for (int i = 0; i < size;i++)
        {
            for (int j = 0;j < sizey;j++)
            {
                for (int k = 0; k < sizez;k++)
                {
                    matriz[i*sizez*sizey+ j*sizez+ k] = 1;
                    matriz1[i*sizez*sizey+ j*sizez+ k] = 1;
                }
            }
        }
        
        ComputeBuffer computeBuffer = new ComputeBuffer(matriz.Length, sizeof(int));
        computeBuffer.SetData(matriz);
        ComputeBuffer computeBuffer1 = new ComputeBuffer(matriz1.Length, sizeof(int));
        computeBuffer1.SetData(matriz1);
        
        computeShader.SetInt("size",size);
        computeShader.SetInt("sizey",sizey);
        computeShader.SetInt("sizez",sizez);
        computeShader.SetBuffer(0,"m0",computeBuffer);
        computeShader.SetBuffer(0,"m1",computeBuffer1);
        computeShader.Dispatch(0,Math.Max( size/5,1),Math.Max(sizey/5,1),Math.Max(sizez/5,1));
        
        computeBuffer1.GetData(matriz1);
        for (int i = 0; i < size;i++)
        {
            String results = "";
            for (int j = 0;j < sizey;j++)
            {
                for (int k = 0; k < sizez;k++)
                {
                    results += " " + matriz1[i * sizez * sizey + j * sizez + k];

                }
                results += "\n";
            }
            Debug.Log(results);
           
        }


    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
