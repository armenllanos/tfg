using UnityEngine;

public class CubeField : MonoBehaviour
{
    public int numCubesX = 10;
    public int numCubesZ = 10;
    public float spacing = 1f;
    public GameObject cubePrefab;
    
    void Start()
    {
        for (int x = 0; x < numCubesX; x++)
        {
            for (int z = 0; z < numCubesZ; z++)
            {
                Vector3 position = new Vector3(x * spacing, 0f, z * spacing);
                Instantiate(cubePrefab, position, Quaternion.identity, transform);
            }
        }
    }
}


