using UnityEngine;
using System.Collections.Generic;

public class AerodynamicSimulation : MonoBehaviour
{
    List<Particle> particles = new List<Particle>();
    public GameObject sphere;
    public int particlesCount = 50; 
    public float particleRadius = 0.2f;
    public float particleMass = 0.00125f;
    public float spawnRadius = 5f;
    public Vector3 gravityDirection = Vector3.down;
    public float windSpeed = 10f;
    
    private float timeCounter = 0f;
    public float lift = 1.0f;
    public float drag = 0.47f; 
    public float airDensity = 1.225f;  
    public Vector3 volumeSize = new Vector3(20, 20, 20);
    private MeshSplitter meshSplitter;
    private List<Triangle> triangles;
    public float springConstant = 50.0f;
    public float dampingConstant = 5.0f;


    void Start()
    {
       
        InitializeParticles();
        meshSplitter = GetComponent<MeshSplitter>();
        if (meshSplitter != null)
        {
            triangles = MeshSplitter.triangles;
        }
        else
        {
            Debug.Log("MeshSplitter component not found.");
        }
    }

    void Update()
    {
        
        timeCounter += Time.deltaTime;
        //InitializeParticles();
        foreach (Particle particle in particles)
        {
            CalculateForces(particle);
            ApplyGravity(particle);
            //ApplySpring(particle);
           
            
        }
        foreach (var particle in particles)
        {
            HandleBoundaryCollisions(particle);
        }
        
        for (int i = 0; i < particles.Count; i++)
    {
        for (int j = i + 1; j < particles.Count; j++)
        {
            ResolveCollisionForParticles(particles[i], particles[j]);
        }

        if (triangles != null)
        {
            ResolveCollisionForParticlesAndTriangles(particles[i], triangles);
        }
    }
    }
    void InitializeParticles()
    {
        int maxAttempts = 1000; 
        

        for (int i = 0; i < particlesCount; i++)
        {
            Vector3 randomVector;
            bool validPosition;

            int attempts = 0;

            do
            {
                validPosition = true;
                randomVector = new Vector3(
                    UnityEngine.Random.Range(-spawnRadius, spawnRadius),
                    UnityEngine.Random.Range(-spawnRadius, spawnRadius),
                    UnityEngine.Random.Range(-spawnRadius, spawnRadius)
                );

                foreach (Particle existingParticle in particles)
                {
                    if (Vector3.Distance(randomVector, existingParticle.position) < 2 * particleRadius)
                    {
                        validPosition = false;
                        break;
                    }
                }

                attempts++;
                if (attempts >= maxAttempts)
                {
                    Debug.LogError("Could not find a valid position for new particle.");
                    return;
                }
            } while (!validPosition);

         

            GameObject oneParticle = Instantiate(sphere, randomVector, Quaternion.identity);
            Particle newParticle = oneParticle.AddComponent<Particle>();
            newParticle.position = randomVector;
            newParticle.velocity = Vector3.zero;
            newParticle.sphere = oneParticle;
            particles.Add(newParticle);
        }
    }
    void CalculateForces(Particle particle)
    {
        Vector3 relativeVelocity = particle.velocity - new Vector3(windSpeed, 0, 0);
        float speed = relativeVelocity.magnitude;
        if (speed == 0)
        {
            return;
        }
        Vector3 dragForce = -0.5f * airDensity * speed * speed * drag * Mathf.PI * Mathf.Pow(particleRadius, 2) * relativeVelocity.normalized;

        
        Vector3 liftDirection = Vector3.Cross(relativeVelocity.normalized, Vector3.up).normalized;
        if (liftDirection == Vector3.zero)
        {
            return; 
        }
        Vector3 liftForce = 0.5f * airDensity * speed * speed * lift * Mathf.PI * Mathf.Pow(particleRadius, 2) * liftDirection;
        
        Vector3 totalForce = dragForce + liftForce;
        
        float damping = 0.99f;
        particle.velocity = particle.velocity * damping + totalForce / particleMass * Time.deltaTime;

        

        particle.position += particle.velocity * Time.deltaTime;
        particle.sphere.transform.position = particle.position;
        
    }

    void ApplyGravity(Particle particle)
    {
        float scaledGravity = 0.1f;
        particle.velocity += gravityDirection * scaledGravity * Time.deltaTime;
        particle.position += particle.velocity * Time.deltaTime;
        particle.sphere.transform.position = particle.position;
    }

    
    void ResolveCollisionForParticles(Particle p1, Particle p2)
{
    Vector3 delta = p1.position - p2.position;
    float distance = delta.magnitude;
    float minDistance = 2 * particleRadius;

    if (distance < minDistance)
    {
        Vector3 normal = delta.normalized;
        float penetrationDepth = minDistance - distance;

      
        Vector3 separation = normal * (penetrationDepth / 2);
        p1.position += separation;
        p2.position -= separation;

       
        p1.sphere.transform.position = p1.position;
        p2.sphere.transform.position = p2.position;

        
        Vector3 relativeVelocity = p1.velocity - p2.velocity;
        float separatingVelocity = Vector3.Dot(relativeVelocity, normal);
        if (separatingVelocity < 0)
        {
            float totalMass = particleMass * 2;
            float impulse = -(1.0f + 0.5f) * separatingVelocity / totalMass;

            Vector3 impulseVector = impulse * normal;
            p1.velocity += impulseVector * (particleMass / totalMass);
            p2.velocity -= impulseVector * (particleMass / totalMass);
        }
    }
}

    void HandleBoundaryCollisions(Particle particle)
    {
        Vector3 minBounds = -volumeSize / 2;
        Vector3 maxBounds = volumeSize / 2;
        float dampingFactor = 0.8f;

        for (int i = 0; i < 3; i++)
        {
            if (particle.position[i] < minBounds[i])
            {
                particle.position[i] = minBounds[i] + 0.01f;
                particle.velocity[i] *= -dampingFactor;
            }
            else if (particle.position[i] > maxBounds[i])
            {
                particle.position[i] = maxBounds[i] - 0.01f;
                particle.velocity[i] *= -dampingFactor;
            }
        }
        particle.transform.position = particle.position;
        
    }
    void ResolveCollisionForParticlesAndTriangles(Particle particle, List<Triangle> triangles)
{
    foreach (var triangle in triangles)
    {
        Vector3 closestPoint = ClosestPointOnTriangle(particle.position, triangle);
        Vector3 delta = particle.position - closestPoint;
        float distance = delta.magnitude;

        if (distance < particleRadius)
        {
            Vector3 normal = delta.normalized;
            float penetrationDepth = particleRadius - distance;

            // Adjust particle position to resolve penetration
            particle.position += normal * penetrationDepth;
            particle.sphere.transform.position = particle.position;

            // Calculate the relative velocity
            Vector3 relativeVelocity = particle.velocity;
            float separatingVelocity = Vector3.Dot(relativeVelocity, normal);

            if (separatingVelocity < 0)
            {
                float impulse = -(1.0f + 0.5f) * separatingVelocity;
                Vector3 impulseVector = impulse * normal;

                // Apply impulse to the particle's velocity
                particle.velocity += impulseVector;
            }
        }
    }
}

Vector3 ClosestPointOnTriangle(Vector3 point, Triangle triangle)
{
    Vector3[] vertices = triangle.GetVertices();
    Vector3 v0 = vertices[1] - vertices[0];
    Vector3 v1 = vertices[2] - vertices[0];
    Vector3 v2 = point - vertices[0];

    float d00 = Vector3.Dot(v0, v0);
    float d01 = Vector3.Dot(v0, v1);
    float d11 = Vector3.Dot(v1, v1);
    float d20 = Vector3.Dot(v2, v0);
    float d21 = Vector3.Dot(v2, v1);

    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    u = Mathf.Clamp01(u);
    v = Mathf.Clamp01(v);
    w = Mathf.Clamp01(w);

    return vertices[0] * u + vertices[1] * v + vertices[2] * w;
}
    void ApplySpring(Particle particle)
    {
       
        Vector3 targetPosition = particle.position + new Vector3(
            Random.Range(-0.05f, 0.05f),
            Random.Range(-0.05f, 0.05f),
            Random.Range(-0.05f, 0.05f)
        );


        Vector3 displacement = targetPosition - particle.position; 
        Vector3 springForce = springConstant * displacement; 
        Vector3 dampingForce = -dampingConstant * particle.velocity; 

        Vector3 force = springForce + dampingForce;

        particle.velocity += force / particleMass * Time.deltaTime; 
        particle.position += particle.velocity * Time.deltaTime; 
        particle.sphere.transform.position = particle.position; 
    }   



    
}