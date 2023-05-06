using System;
using System.Collections;
using System.Collections.Generic;
using OpenCover.Framework.Model;
using UnityEngine;
using UnityEngine.UIElements;

public class pruebaInterfaz : MonoBehaviour
{
    public UIDocument menu;

    private void OnEnable()
    {

        menu = GetComponent<UIDocument>();
        VisualElement root = menu.rootVisualElement;
        Label texto = root.Q<Label>(className:"clase");
        Label texto2 = root.Q<Label>("nombre elemento");
        
    }

    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
