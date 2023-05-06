using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using DefaultNamespace;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Accessibility;
using UnityEngine.UIElements;
using Debug = UnityEngine.Debug;
using Random = UnityEngine.Random;


public class Creator2D : MonoBehaviour
{
    public ComputeShader computeShader;

    public RenderTexture renderTexture;

    // Start is called before the first frame update
    public int numCubesX = 10;
    public int numCubesY = 10;
    public float spacing = 1f;
    public GameObject cubePrefab;
    public bool showPressure = false;
    public float speed1;
    public float speed2;
    public int numberOfFans;
    public float cellSize;
    public Camera myCamera;
    public GameObject fanModel;

    private Dictionary<int, double> thrusts = new Dictionary<int, double>();
    private Dictionary<int, double> powers = new Dictionary<int, double>();
    private GameObject[,] cubes;
    private Fluid2D fluid;
    private ArrayList gameObjects = new ArrayList();
    private double[] speeds;
    private int j = 0;
    public bool GPU = false;
    float opacity = 0;
    UIDocument menu;
    private double[,] smoke;
    private bool started;
    private int axis;
    private int fansX;
    private int fansY;
    private bool selection;
    double totalThrust = 0, totalCurrent = 0;
    Stopwatch stopwatch;

    private Button startButton;
    private SliderInt numFansX;
    private SliderInt numFansY;
    private Label labelNumFansX;
    private Label labelNumFansY;
    private Label frameLabel;

    private Ventilador[,] fans;
    double currentX = -0.09194635; //Data calculated by the script
    double currentX2 = 0.03139418;
    double currentCutY = 0.19078497667231487;

    private const double voltaje = 12;
    private const double MAX_VELOCITY = 9.2;

    double thrustX = -1.34674003e-06;
    double thrustX2 = 4.02562367e-09;
    double thrustCutY = 0.002088167064807224;

    double RPMSlope = 479.02721321;
    double RPMCutY = 438.0384358478027;


    private void OnEnable()
    {
        menu = GetComponent<UIDocument>();
        VisualElement root = menu.rootVisualElement;
        root.Q<VisualElement>("centro").visible = false;

        numFansX = root.Q<SliderInt>("SliderNumberFansX");
        numFansY = root.Q<SliderInt>("SliderNumberFansY");

        frameLabel = root.Q<Label>("fps");

        labelNumFansX = root.Q<Label>("numberFansX");
        labelNumFansY = root.Q<Label>("numberFansY");

        startButton = root.Q<Button>("startButton");
        startButton.text = "Next";
        IStyle color = startButton.style;
        color.backgroundColor = new StyleColor(new Color((float)0.37,(float)0.78,(float)0.94));
        numFansX.RegisterCallback<ClickEvent, string>(updateNumberOfFans, "X");
        numFansY.RegisterCallback<ClickEvent, string>(updateNumberOfFans, "Y");
        startButton.RegisterCallback<ClickEvent>(chooseSize);
        fans = new Ventilador[16,16];
        root.Q<VisualElement>("datos").Q<VisualElement>("totalThrust").Q<Label>("thrust").text = 0.ToString();
        root.Q<VisualElement>("datos").Q<VisualElement>("totalPower").Q<Label>("power").text = 0.ToString();
        
        root.Q<Button>("configureButton").RegisterCallback<ClickEvent>(showConfiguration);
        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 16; j++)
            {
                fans[i,j] = new Ventilador(i + 1, 0, calclateRPM(0), calculateThrust(0.0), calculateCurrent(0.0));
                root.Q<VisualElement>("fanColumn" + (i + 1)).Q<Button>("fan"+(j+1)).RegisterCallback<ClickEvent,Ventilador>(showFan,fans[i,j]);
                root.Q<VisualElement>("fanColumn" + (i + 1)).Q<Button>("fan" + (j + 1)).visible = false;
                
            }
        }
    }

    private void showConfiguration(ClickEvent evt)
    {
        VisualElement root = menu.rootVisualElement;
        for (int i = 0; i < fansX; i++)
        {
            for (int j = 0; j < fansY; j++)
            {
                root.Q<VisualElement>("fanColumn" + (i + 1)).Q<Button>("fan" + (j + 1)).visible = true;
            }
        }
        
    }

    private void showFan(ClickEvent evt, Ventilador fan)
    {
        VisualElement root = menu.rootVisualElement;
        int row, column;
        row = Math.Abs( position.Item2 - 16);
        column = Math.Abs(position.Item1 - 16);
        root.Q<VisualElement>("fanInfo").Q<Label>("number").text =
            "Fan of row: " + position.Item2 + " column: " + position.Item1;
        root.Q<VisualElement>("fanInfo").Q<Label>("Velocity").text =  fans[column,row].actualVelocity.ToString()+"m/s";
        root.Q<VisualElement>("fanInfo").Q<Label>("RPM").text = fans[column,row].actualRPM.ToString().Substring(0, Math.Min(6, fans[column,row].actualRPM.ToString().Length));
        root.Q<VisualElement>("fanInfo").Q<Label>("Amperaje").text = fans[column,row].actualCurrent.ToString().Substring(0, Math.Min(6, fans[column,row].actualCurrent.ToString().Length))+"A";
        root.Q<VisualElement>("fanInfo").Q<Label>("Empuje").text = fans[column,row].actualThrust.ToString().Substring(0, Math.Min(6, fans[column,row].actualThrust.ToString().Length))+"N";
        root.Q<VisualElement>("fanInfo").Q<SliderInt>("VelocityPercentaje").value =
            (int)(fans[column,row].actualVelocity / MAX_VELOCITY)*100 ;
        root.Q<VisualElement>("fanInfo").visible = true;
        root.Q<SliderInt>("VelocityPercentaje")
            .RegisterCallback<ClickEvent, Ventilador>(changeVelocity, fans[i,j]);
    }
    private void changeVelocity(ClickEvent evt, Ventilador fan)
    {
        VisualElement root = menu.rootVisualElement;
        int row, column;
        row = Math.Abs( number.Item2 - 16);
        column = Math.Abs(number.Item1 - 16);
        double velocity, RPM, current, thrust;
        double oldThrust = fans[column,row].actualThrust;
        double oldCurrent = fans[column,row].actualCurrent;
        int percentaje = root.Q<SliderInt>("VelocityPercentaje").value;
        velocity = percentaje * MAX_VELOCITY / 100;
        RPM = calclateRPM(velocity);
        current = calculateCurrent(velocity);
        thrust = calculateThrust(RPM);
        fans[column,row].actualThrust = thrust;
        fans[column,row].actualCurrent = current;
      
        String positionText = root.Q<VisualElement>("fanInfo").Q<Label>("number").text;
        String[] tokens = positionText.Split(" ");
       


        root.Q<VisualElement>("fanInfo").Q<Label>("Velocity").text =  (velocity).ToString()+"m/s";
        fans[column, row].actualVelocity = velocity;
        
        root.Q<VisualElement>("fanInfo").Q<Label>("RPM").text = RPM.ToString().Substring(0, Math.Min(6, RPM.ToString().Length));
        fans[column, row].actualRPM = RPM;
        
        root.Q<VisualElement>("fanInfo").Q<Label>("Amperaje").text = current.ToString().Substring(0, Math.Min(6, current.ToString().Length))+"A";
        fans[column, row].actualCurrent = current;
        
        root.Q<VisualElement>("fanInfo").Q<Label>("Empuje").text = thrust.ToString().Substring(0, Math.Min(6, thrust.ToString().Length))+"N";
        fans[column, row].actualThrust = thrust;
        if (!selection)
        {
            if (column == axis)
            {
                if (fluid != null)
                {
                    fluid.speed[row] = velocity;
                }
                this.speeds[row] = velocity;
            }
        }


        totalCurrent -= oldCurrent;
        totalCurrent += current;

        totalThrust -= oldThrust;
        totalThrust += thrust;

        root.Q<VisualElement>("totalThrust").Q<Label>("thrust").text =
            (totalThrust).ToString().Substring(0, Math.Min(6, thrust.ToString().Length));
        root.Q<VisualElement>("totalPower").Q<Label>("power").text =
            (totalCurrent).ToString().Substring(0, Math.Min(6, thrust.ToString().Length));
    }

    private void updateNumberOfFans(ClickEvent evt, string pos)
    {
        if (pos.Equals("X"))
        {
            labelNumFansX.text = numFansX.value.ToString();
        }
        else
        {
            labelNumFansY.text = numFansY.value.ToString();
        }
    }

    private void chooseSize(ClickEvent evt)
    {
        numFansY.visible = false;
        numFansY.focusable = false;
        numFansX.highValue = numFansX.value;
        numberOfFans = numFansY.value;
        fansX = numFansX.value;
        fansY = numFansY.value;
        labelNumFansY.visible = false;
        numCubesY *= numFansY.value;
        myCamera.transform.SetPositionAndRotation(
            new Vector3(numCubesX / 2, numCubesY / 2, -Math.Max(numCubesX, numCubesY)), Quaternion.identity);
        cubes = new GameObject[numCubesX, numCubesY];
        for (int x = 0; x < numCubesX; x++)
        {
            for (int y = 0; y < numCubesY; y++)
            {
                Vector3 position = new Vector3(x * spacing + (float)0.5, y * spacing + (float)0.5,
                    0 * spacing + (float)0.5);
                cubes[x, y] = Instantiate(cubePrefab, position, Quaternion.identity, transform);
                gameObjects.Add(cubes[x, y]);
            }
        }

        VisualElement root = menu.rootVisualElement;
        root.Q<Label>("Ventiladores").text = "Elige el eje de corte";
        root.Q<Label>("Y").visible = false;
        startButton.UnregisterCallback<ClickEvent>(chooseSize);
        startButton.RegisterCallback<ClickEvent>(configure);
    }

    private void configure(ClickEvent evt)
    {
        
        VisualElement root = menu.rootVisualElement;
        axis = numFansX.value;
        for (int i = 0; i < fansX; i++)
        {
            for (int j = 0; j < fansY; j++)
            {
                root.Q<VisualElement>("fanColumn" + (i + 1)).Q<Button>("fan" + (j + 1)).visible = true;
            }
        }

        startButton.UnregisterCallback<ClickEvent>(configure);
        startButton.RegisterCallback<ClickEvent>(start);
    }

    private void start(ClickEvent evt)
    {
        if (fluid == null)
        {
            fluid = new Fluid2D(1000, numCubesX, numCubesY, cellSize, speeds, "x0", computeShader,
                numCubesY / (numberOfFans + 1), numberOfFans);
        }

        for (int i = 0; i < numberOfFans; i++)
        {
            Vector3 position = new Vector3(-5, i * 20 + 10,
                10);
            Vector3 rotation = new Vector3(0, -20+(numFansY.value*2), 0);
            GameObject fan = Instantiate(fanModel, position, Quaternion.identity, transform);
            fan.transform.Rotate(rotation);
            gameObjects.Add(fan);
        }

        numFansX.visible = false;
        numFansX.focusable = false;

        numFansY.visible = false;
        numFansY.focusable = false;
        started = !started;
        if (startButton.text.Equals("Start"))
        {
            startButton.text = "Stop";
        }
        else
        {
            startButton.text = "Start";
        }
    }

    private double calculateCurrent(double d)
    {
        if (d == 0)
        {
            return 0;
        }

        return currentX2 * Math.Pow(d, 2) + d * currentX + currentCutY;
    }

    private double calculateThrust(double d)
    {
        if (d == 0)
        {
            return 0;
        }

        return thrustX2 * Math.Pow(d, 2) + d * thrustX + thrustCutY;
    }

    private double calclateRPM(double i)
    {
        if (i == 0)
        {
            return 0;
        }

        return RPMSlope * i + RPMCutY;
    }

    void Start()
    {
        started = false;
        myCamera.transform.SetPositionAndRotation(new Vector3(numCubesX / 2, numCubesY / 2, -numCubesX),
            Quaternion.identity);
        speeds = new double[8];
    }

    // Update is called once per frame
    void Update()
    {
        if (started)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            Color myColor;
            if (GPU)
            {
                fluid.simulateGPU(Time.deltaTime);
            }
            else
            {
                fluid.simulate(Time.deltaTime);
            }

            double minPressure = fluid.p[0, 0];
            double maxPressure = fluid.p[0, 0];
            for (int q = 0; q < numCubesX + 2; q++)
            {
                for (int w = 0; w < numCubesY + 2; w++)
                {
                    minPressure = Math.Min(minPressure, fluid.p[q, w]);
                    maxPressure = Math.Max(maxPressure, fluid.p[q, w]);
                }
            }

            for (int x = 0; x < numCubesX; x++)
            {
                for (int y = 0; y < numCubesY; y++)
                {
                    if (showPressure)
                    {
                        Color finalColor = getColorRGB(fluid.p[x + 1, y + 1], minPressure, maxPressure);
                        cubes[x, y].GetComponent<Renderer>().material.color =
                            finalColor;
                    }
                    else
                    {
                        
                        Color finalColor = new Color((float)fluid.smokeField[x + 1, y + 1] / (float)Fluid2D.MAX_SPEED,
                            (float)fluid.smokeField[x + 1, y + 1] / (float)Fluid2D.MAX_SPEED,
                            (float)fluid.smokeField[x + 1, y + 1] / (float)Fluid2D.MAX_SPEED,
                            (float)Math.Abs(fluid.smokeField[x + 1, y + 1] - MAX_VELOCITY) /
                            (float)Fluid2D.MAX_SPEED); //getColorRGB(fluid.smokeField[x + 1, y + 1],0,1);
                        /*Debug.Log("opacity: " + fluid.smokeField[x, y, z]);*/
                        cubes[x, y].GetComponent<Renderer>().material.color =
                            finalColor;
                        /*Debug.Log("opacity: " + fluid.smokeField[x, y, z] / speed + "smoke value" +
                                 fluid.smokeField[x, y, z] / speed);*/
                        if (fluid.s[x, y] == 0)
                        {
                            
                        }
                    }
                }
            }

            stopwatch.Stop();
            long wat = stopwatch.ElapsedMilliseconds;
            double time = 1.0 / (stopwatch.ElapsedMilliseconds/1000.0);
            frameLabel.text = time.ToString().Substring(0, Math.Min(3, time.ToString().Length));
        }
    }

    private Color getColorRGB(double p, double min, double max)
    {
        double val = Math.Min(Math.Max(p, min), max - 0.0001);
        double d = max - min;
        val = d == 0.0 ? 0.5 : (val - min) / d;
        double m = 0.25;
        double num = Math.Floor(val / m);
        double s = (val - num * m) / m;
        double r = 0, g = 0, b = 0;

        switch (num)
        {
            case 0:
                r = 0.0;
                g = s;
                b = 1.0;
                break;
            case 1:
                r = 0.0;
                g = 1.0;
                b = 1.0 - s;
                break;
            case 2:
                r = s;
                g = 1.0;
                b = 0.0;
                break;
            case 3:
                r = 1.0;
                g = 1.0 - s;
                b = 0.0;
                break;
        }

        double mean = r + g + b;
        mean /= 3;
        return new Color((float)r, (float)g, (float)b, 1);
    }

    private Color getColor(double val)
    {
        val = Math.Min(Math.Max(val, 0), 1 - 0.0001);
        int d = 1;
        //val = d == 0.0 ? 0.5 : (val - 0) / d;
        double m = 0.25;
        int num = (int)Math.Floor(val / m);
        double s = (val - num * m) / m;

        return new Color((float)0, (float)0, (float)0, (float)((1 - val)));
    }


    public static void LinearRegression(double[] xVals, double[] yVals,
        int inclusiveStart, int exclusiveEnd,
        out double rsquared, out double yintercept,
        out double slope)
    {
        Debug.Assert(xVals.Length == yVals.Length);
        double sumOfX = 0;
        double sumOfY = 0;
        double sumOfXSq = 0;
        double sumOfYSq = 0;
        double ssX = 0;
        double ssY = 0;
        double sumCodeviates = 0;
        double sCo = 0;
        double count = exclusiveEnd - inclusiveStart;

        for (int ctr = inclusiveStart; ctr < exclusiveEnd; ctr++)
        {
            double x = xVals[ctr];
            double y = yVals[ctr];
            sumCodeviates += x * y;
            sumOfX += x;
            sumOfY += y;
            sumOfXSq += x * x;
            sumOfYSq += y * y;
        }

        ssX = sumOfXSq - ((sumOfX * sumOfX) / count);
        ssY = sumOfYSq - ((sumOfY * sumOfY) / count);
        double RNumerator = (count * sumCodeviates) - (sumOfX * sumOfY);
        double RDenom = (count * sumOfXSq - (sumOfX * sumOfX))
                        * (count * sumOfYSq - (sumOfY * sumOfY));
        sCo = sumCodeviates - ((sumOfX * sumOfY) / count);

        double meanX = sumOfX / count;
        double meanY = sumOfY / count;
        double dblR = RNumerator / Math.Sqrt(RDenom);
        rsquared = dblR * dblR;
        yintercept = meanY - ((sCo / ssX) * meanX);
        slope = sCo / ssX;
    }
}