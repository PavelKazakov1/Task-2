using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;

class Program
{
    static void Main(string[] args)
    {
        // Constants for signal generation
        const bool muteNoise = false;
        const int Fs = 8192; 
        const double T = 1;
        const int N = Fs * (int)T;
        var t = Enumerable.Range(0, N).Select(i => (double)i / Fs).ToArray(); 

        // Signal parameters
        const double f1 = 40, f2 = 100, A1 = 1, A2 = 0.5;

        // Generate the signal
        var signal = t.Select(ti => A1 * Math.Sin(2 * Math.PI * f1 * ti) + A2 * Math.Sin(2 * Math.PI * f2 * ti)).ToArray();

        // Parameters for noise generation
        var mean = 0.0;
        var variance = 1.0;

        // Generate noise
        var noise = t.Select(_ => (!muteNoise ? Math.Sqrt(-2.0 * Math.Log(Random.NextDouble())) * Math.Cos(2.0 * Math.PI * Random.NextDouble()) * Math.Sqrt(variance) + mean : 0.0)).ToArray();

        // Combine signal and noise
        var signalPlusNoise = signal.Zip(noise, (s, n) => s + n).ToArray();

        // Size of the moving average window
        const int windowSize = 5;

        // Apply moving average filter to signal, noise, and signal+noise
        var filteredSignal = MovingAverageFilter(signal, windowSize);
        var filteredNoise = MovingAverageFilter(noise, windowSize);
        var filteredSignalPlusNoise = MovingAverageFilter(signalPlusNoise, windowSize);

        // Calculate SNR  before and after filtering
        var snrBeforeFiltering = CalculateSNR(signal, noise);
        var snrAfterFiltering = CalculateSNR(filteredSignal, filteredNoise);

        // Print SNR values
        Console.WriteLine($"SNR before filtering: {snrBeforeFiltering}, SNR after filtering: {snrAfterFiltering}");

        // Perform FFT 
        CalculateFFT(signal, Fs);
        CalculateFFT(noise, Fs);
        CalculateFFT(signalPlusNoise, Fs);
        CalculateFFT(filteredSignal, Fs);
        CalculateFFT(filteredNoise, Fs);
        CalculateFFT(filteredSignalPlusNoise, Fs);

        // Calculate autocorrelation of signal, noise, and signal+noise
        var acfSignal = Autocorrelation(signal);
        var acfNoise = Autocorrelation(noise);
        var acfSignalPlusNoise = Autocorrelation(signalPlusNoise);

        Console.WriteLine("Results:");
        Console.WriteLine($"SNR before filtering: {snrBeforeFiltering}, SNR after filtering: {snrAfterFiltering}");
        Console.WriteLine("Calculation Completed.");
        Console.ReadLine();
    }

    static Random Random = new Random();

    // method to perform FFT
    static void CalculateFFT(double[] data, int Fs)
    {
        // Convert data to Complex32 array
        var complexData = data.Select(d => new Complex32((float)d, 0)).ToArray();

        // Perform FFT on the complex data
        Fourier.Forward(complexData);
    }

    // method to apply moving average filter
    static double[] MovingAverageFilter(double[] data, int windowSize)
    {
        var kernel = Enumerable.Repeat(1.0 / windowSize, windowSize).ToArray();
        var filteredData = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
        {
            for (int k = 0; k < windowSize; k++)
            {
                filteredData[i] += (i - k >= 0 ? data[i - k] : 0) * kernel[k];
            }
        }
        return filteredData;
    }

    // method to calculate autocorrelation
    static double[] Autocorrelation(double[] data)
    {
        var result = new double[data.Length];
        for (int lag = 0; lag < data.Length; lag++)
        {
            for (int i = 0; i < data.Length - lag; i++)
            {
                result[lag] += data[i] * data[i + lag];
            }
            result[lag] /= data.Length;
        }
        return result.Take(data.Length / 2).ToArray();
    }

    // method to calculate SNR (Signal-to-Noise Ratio)
    static double CalculateSNR(double[] signal, double[] noise)
    {
        var powerSignal = signal.Select(val => val * val).Sum() / signal.Length;
        var powerNoise = noise.Select(val => val * val).Sum() / noise.Length;
        return 10 * Math.Log10(powerSignal / powerNoise);
    }
}
