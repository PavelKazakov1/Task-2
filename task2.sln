
using System;

public class EuclideanDistance
{

    public static double Calculate(double x1, double y1, double x2, double y2)
    {
        return Math.Sqrt(Math.Pow(x2 - x1, 2) + Math.Pow(y2 - y1, 2));
    }

    public static double Loss(double x, double y, double[][] anchors, double[] distances)
    {
        double loss = 0;
        for (int i = 0; i < anchors.Length; i++)
        {
            double measuredDistance = distances[i];
            double calculatedDistance = Calculate(x, y, anchors[i][0], anchors[i][1]);
            loss += Math.Pow(calculatedDistance - measuredDistance, 2);
        }
        return loss;
    }

    public static double[] Gradient(double x, double y, double[][] anchors, double[] distances)
    {
        double gradX = 0, gradY = 0;
        for (int i = 0; i < anchors.Length; i++)
        {
            double dx = x - anchors[i][0];
            double dy = y - anchors[i][1];
            double measuredDistance = distances[i];
            double calculatedDistance = Calculate(x, y, anchors[i][0], anchors[i][1]);
            gradX += (calculatedDistance - measuredDistance) * (dx / calculatedDistance);
            gradY += (calculatedDistance - measuredDistance) * (dy / calculatedDistance);
        }
        return new double[] { gradX, gradY };
    }

    public static double[] GradientDescent(double[][] anchors, double[] distances, double[] initialGuess, double learningRate, int iterations)
    {
        double x = initialGuess[0];
        double y = initialGuess[1];
        for (int i = 0; i < iterations; i++)
        {
            double[] gradients = Gradient(x, y, anchors, distances);
            x -= learningRate * gradients[0];
            y -= learningRate * gradients[1];
        }
        return new double[] { x, y };
    }

    public static void Main(string[] args)
    {
        double[][] anchors = new double[][] { new double[] { 1, 1 }, new double[] { 2, 3 }, new double[] { 4, 2 } };
        double[] distances = new double[] { 2.236, 2.828, 2.236 };
        double[] initialGuess = new double[] { 0, 0 };
        double learningRate = 0.01;
        int iterations = 1000;
        double[] result = GradientDescent(anchors, distances, initialGuess, learningRate, iterations);
        Console.WriteLine($"Estimated coordinates: ({result[0]}, {result[1]})");
    }
}
