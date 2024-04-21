using System;
using System.Collections.Generic;
using System.Linq;

class Program
{
    static void Main(string[] args)
    {
        // Constants
        double[,] anchors = { { 0, 0 }, { 0, 10 }, { 10, 0 } };
        double[] objectCoordinates = { 5, 5 };
        double noiseVariance = 1;
        int iterations = 1000;
        double learningRate = 0.01;

        // Distances to base stations
        double[] distances = new double[anchors.GetLength(0)];
        for (int i = 0; i < anchors.GetLength(0); i++)
        {
            distances[i] = EuclideanDistance(objectCoordinates[0], objectCoordinates[1], anchors[i, 0], anchors[i, 1]);
        }

        // Angles for triangulation (not used in C# version)
        double[] angles = new double[anchors.GetLength(0)];
        for (int i = 0; i < anchors.GetLength(0); i++)
        {
            angles[i] = CalculateAngle(objectCoordinates[0], objectCoordinates[1], anchors[i, 0], anchors[i, 1]);
        }

        // Trilateration without noise
        double[] trilaterationResult = Trilateration(anchors, distances);

        // Add noise to distances
        double[] noisyDistances = AddNoise(distances, noiseVariance);

        // Trilateration with noise
        double[] noisyTrilaterationResult = Trilateration(anchors, noisyDistances);

        // Gradient descent
        double[] gradientDescentResult = GradientDescent(anchors, distances, iterations, learningRate);

        // Print results
        Console.WriteLine("Initial Data:");
        Console.WriteLine("Base Stations Coordinates:");
        PrintArray(anchors);
        Console.WriteLine("Object Coordinates:");
        PrintArray(objectCoordinates);

        Console.WriteLine("\nCalculations:");
        Console.WriteLine("Distances to Base Stations:");
        PrintArray(distances);
        Console.WriteLine("Trilateration Result (without noise):");
        PrintArray(trilaterationResult);

        Console.WriteLine("\nExperiment with Noise:");
        Console.WriteLine("Noisy Distances:");
        PrintArray(noisyDistances);
        Console.WriteLine("Trilateration Result (with noise):");
        PrintArray(noisyTrilaterationResult);

        Console.WriteLine("\nGradient Descent Result:");
        PrintArray(gradientDescentResult);
    }

    // Function to calculate Euclidean distance
    static double EuclideanDistance(double x1, double y1, double x2, double y2)
    {
        return Math.Sqrt(Math.Pow(x2 - x1, 2) + Math.Pow(y2 - y1, 2));
    }

    // Function to calculate angle
    static double CalculateAngle(double x1, double y1, double x2, double y2)
    {
        return Math.Atan2(y2 - y1, x2 - x1);
    }

    // Function for trilateration
    static double[] Trilateration(double[,] anchors, double[] distances)
    {
        double x1 = anchors[0, 0];
        double y1 = anchors[0, 1];
        double x2 = anchors[1, 0];
        double y2 = anchors[1, 1];
        double x3 = anchors[2, 0];
        double y3 = anchors[2, 1];

        double d1 = distances[0];
        double d2 = distances[1];
        double d3 = distances[2];

        double A = 2 * (x2 - x1);
        double B = 2 * (y2 - y1);
        double C = 2 * (x3 - x1);
        double D = 2 * (y3 - y1);

        double E = Math.Pow(d1, 2) - Math.Pow(d2, 2) - Math.Pow(x1, 2) + Math.Pow(x2, 2) - Math.Pow(y1, 2) + Math.Pow(y2, 2);
        double F = Math.Pow(d1, 2) - Math.Pow(d3, 2) - Math.Pow(x1, 2) + Math.Pow(x3, 2) - Math.Pow(y1, 2) + Math.Pow(y3, 2);

        double x = (E * D - F * B) / (A * D - B * C);
        double y = (F * A - E * C) / (A * D - B * C);

        return new double[] { x, y };
    }

    // Function to add noise to distances
    static double[] AddNoise(double[] distances, double variance)
    {
        Random random = new Random();
        return distances.Select(distance => distance + Math.Sqrt(variance) * (random.NextDouble() - 0.5)).ToArray();
    }

    // Function for gradient descent
    static double[] GradientDescent(double[,] anchors, double[] distances, int iterations, double learningRate)
    {
        double[] initialGuess = { 1, 1 };
        double x = initialGuess[0];
        double y = initialGuess[1];

        for (int i = 0; i < iterations; i++)
        {
            double[] gradient = CalculateGradient(x, y, anchors, distances);
            x -= learningRate * gradient[0];
            y -= learningRate * gradient[1];
        }

        return new double[] { x, y };
    }

    // Function to calculate gradient
    static double[] CalculateGradient(double x, double y, double[,] anchors, double[] distances)
    {
        double gradX = 0, gradY = 0;

        for (int i = 0; i < anchors.GetLength(0); i++)
        {
            double dx = x - anchors[i, 0];
            double dy = y - anchors[i, 1];
            double measuredDistance = distances[i];
            double calculatedDistance = EuclideanDistance(x, y, anchors[i, 0], anchors[i, 1]);

            gradX += (calculatedDistance - measuredDistance) * (dx / calculatedDistance);
            gradY += (calculatedDistance - measuredDistance) * (dy / calculatedDistance);
        }

        return new double[] { gradX, gradY };
    }

    // Function to print array
    static void PrintArray(double[] array)
    {
        Console.WriteLine("[" + string.Join(", ", array.Select(x => x.ToString("0.##"))) + "]");
    }

    // Function to print 2D array
    static void PrintArray(double[,] array)
    {
        for (int i = 0; i < array.GetLength(0); i++)
        {
            PrintArray(Enumerable.Range(0, array.GetLength(1)).Select(j => array[i, j]).ToArray());
        }
    }
}
