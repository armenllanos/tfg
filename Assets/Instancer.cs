using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Instancer : MonoBehaviour
{

    public int Instances;
    public Mesh mesh;
    public Material[] Materials;
    private List<List<Matrix4x4>> Batches = new List<List<Matrix4x4>>();

    private void RenderBatches()
    {
        foreach(List<Matrix4x4> Batch in Batches){
            for(int i = 0; i < mesh.subMeshCount; i++)
            {
                Graphics.DrawMeshInstanced(mesh, i, Materials[i], Batch);

            }
        }

    }

    void Update()
    {
        RenderBatches();
    }

    // Start is called before the first frame update
    void Start()
    {
        int AddedMatricies = 0;
        for(int i = 0; i < Instances; i++)
        {
            if(AddedMatricies < 1000)
            {
                Batches[Batches.Count - 1].Add(Matrix4x4.TRS(new Vector3(Random.Range(0, 50), Random.Range(0, 50), Random.Range(0, 50)),Random.rotation,Vector3.back));
                AddedMatricies += 1;
            }
            else
            {
                Batches.Add(new List<Matrix4x4>());
                AddedMatricies = 0;
            }
        }
    }

    // Update is called once per frame
}
