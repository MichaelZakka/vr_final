using UnityEngine;
using UnityEngine.SceneManagement;

public class MainMenu : MonoBehaviour
{
    public void StartScene()
    {
        Time.timeScale = 1f;
    }
    public void X1()
    {
        Time.timeScale = 1.5f;
    }
    public void X2()
    {
        Time.timeScale = 2f;
    }

    public void ResetScene()
    {
        Scene currentScene = SceneManager.GetActiveScene();
        
        SceneManager.LoadScene(currentScene.name);
    }


    public void PauseScene()
    {
        Time.timeScale = 0f;
    }

}
