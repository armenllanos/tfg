using UnityEngine.UIElements;
using Label = System.Reflection.Emit.Label;

namespace DefaultNamespace
{
    public class Ventilador
    {
        public int fanNum;
        public double actualVelocity;
        public double actualRPM;
        public double actualThrust;
        public double actualCurrent;
        public int fanColumn;
        private const double MAX_VELOCITY = 9.2;

        public Ventilador(int fanNum,int fanColumn, double actualVelocity, double actualRpm, double actualThrust, double actualCurrent)
        {
            this.fanNum = fanNum;
            this.fanColumn = fanColumn;
            this.actualVelocity = actualVelocity;
            actualRPM = actualRpm;
            this.actualThrust = actualThrust;
            this.actualCurrent = actualCurrent;
        }
    }
}