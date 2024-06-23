using UnityEngine;
using System.Collections.Generic;

public class MeshSplitter : MonoBehaviour
{
    public GameObject planeObject;
    private Mesh originalMesh;
    private Plane cuttingPlane;
    private Vector3[] vertices;
    private int[] triangles;

    private List<Triangle> triangleList = new List<Triangle>(); // List to store Triangle objects

    void Start()
    {
        originalMesh = GetComponent<MeshFilter>().mesh;
        cuttingPlane = new Plane(planeObject.transform.up, planeObject.transform.position);
        vertices = originalMesh.vertices;
        triangles = originalMesh.triangles;

        SplitMesh(originalMesh);
    }

    void SplitMesh(Mesh mesh)
    {
        List<Vector3> abovePlaneVertices = new List<Vector3>();
        List<Vector3> belowPlaneVertices = new List<Vector3>();

        for (int i = 0; i < triangles.Length; i += 3)
        {
            Vector3 v0 = transform.TransformPoint(vertices[triangles[i]]);
            Vector3 v1 = transform.TransformPoint(vertices[triangles[i + 1]]);
            Vector3 v2 = transform.TransformPoint(vertices[triangles[i + 2]]);

            bool v0Above = cuttingPlane.GetSide(v0);
            bool v1Above = cuttingPlane.GetSide(v1);
            bool v2Above = cuttingPlane.GetSide(v2);

            if (v0Above || v1Above || v2Above)
            {
                abovePlaneVertices.Add(v0);
                abovePlaneVertices.Add(v1);
                abovePlaneVertices.Add(v2);
                triangleList.Add(new Triangle(v0, v1, v2)); // Add triangle to list
            }
            if (!v0Above || !v1Above || !v2Above)
            {
                belowPlaneVertices.Add(v0);
                belowPlaneVertices.Add(v1);
                belowPlaneVertices.Add(v2);
                triangleList.Add(new Triangle(v0, v1, v2)); // Add triangle to list
            }

            Debug.Log($"Added triangle: {v0}, {v1}, {v2}");
        }

        CreateMeshPart(abovePlaneVertices, "AbovePlaneMesh");
        CreateMeshPart(belowPlaneVertices, "BelowPlaneMesh");
    }

    void CreateMeshPart(List<Vector3> verticesList, string name)
    {
        if (verticesList.Count == 0) return;

        GameObject newObject = new GameObject(name, typeof(MeshFilter), typeof(MeshRenderer), typeof(MeshCollider));
        newObject.transform.position = transform.position;
        newObject.transform.rotation = transform.rotation;

        Mesh newMesh = new Mesh();
        newMesh.vertices = verticesList.ToArray();
        
        List<int> newTriangles = new List<int>();
        for (int i = 0; i < verticesList.Count; i += 3)
        {
            newTriangles.Add(i);
            newTriangles.Add(i + 1);
            newTriangles.Add(i + 2);
        }
        newMesh.triangles = newTriangles.ToArray();

        newMesh.RecalculateNormals();
        newMesh.RecalculateBounds();

        newObject.GetComponent<MeshFilter>().mesh = newMesh;
        newObject.GetComponent<MeshCollider>().sharedMesh = newMesh;
        newObject.GetComponent<MeshCollider>().convex = true;

        Debug.Log($"{name} created with {verticesList.Count / 3} triangles");
    }

    public List<Triangle> GetTriangles()
    {
        return triangleList;
    }

    void OnDrawGizmos()
    {
        if (planeObject != null)
        {
            Gizmos.color = Color.red;
            Gizmos.DrawLine(planeObject.transform.position, planeObject.transform.position + planeObject.transform.up * 5);
            Gizmos.DrawLine(planeObject.transform.position, planeObject.transform.position - planeObject.transform.up * 5);
        }
    }
}


public class Triangle
{
    public Vector3 v0;
    public Vector3 v1;
    public Vector3 v2;

    public Triangle(Vector3 v0, Vector3 v1, Vector3 v2)
    {
        this.v0 = v0;
        this.v1 = v1;
        this.v2 = v2;
    }
}
