//-----------------------------------------------------------------------
// Copyright 2018 Eli Lilly and Company
//
// Licensed under the Apache License, Version 2.0 (the "License");
//
// you may not use this file except in compliance with the License.
//
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//-----------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Science.Chemistry;
using SignalProcessing;
using System.Diagnostics;
using System.Threading.Tasks;
using MassSpectrometry;

public class ChargeDetector
{

    SortedList<double, float> mzList = null;

    SortedList<float, double> intensityList = new SortedList<float, double>();

    //DEBUG:  if set to true clusters that are rejected will have thier failure evidence printed to output
    private const bool verbose = true;

    public ChargeDetector(SortedList<double, float> points, float minIntensity = 0)
    {
        mzList = points;
        intensityList.Clear();

        foreach (var aPoint in mzList)
        {
            if (!intensityList.ContainsKey(aPoint.Value) && aPoint.Value > minIntensity) intensityList.Add(aPoint.Value, aPoint.Key);
        }
    }


    public List<Cluster> DetectChargeStates(double MassTolerance, int MinCharge = 1, int MaxCharge = 40, double massRangeMin = 0, double massRangeMax = double.MaxValue, bool applyNoiseFilter = true)
    {

        // Can be Parallel!!!
        var finalResult = new List<Cluster>(); // Dictionary<int, System.Collections.Generic.SortedList<double, float>>();

        foreach (var chargeState in Enumerable.Range(MinCharge, (MaxCharge - MinCharge) + 1).Reverse())
        {
            var result = DetectCharge(chargeState, MassTolerance, massRangeMin, massRangeMax);


            //TODO: filter out peaks that don't have others from the isotope distribution
            //TODO: filter out distributions that have gaps (or missing peaks) -- open up tolerance a bit to be sure though...

            double zSpace = SignalProcessor.ProtonMass / (double)chargeState;  // changed to accomodate mass defect

            //TODO: take the groupings of peaks and fill in the ones missing from the distribution.  

            var result2 = FindIsotopeGroups(result, chargeState, MassTolerance, applyNoiseFilter);

            foreach (var aCluster in result2)
            {
                bool skip = false;

                int aClusterCorePeakCount = aCluster.Peaks.Count(p => p.IsCorePeak);

                // go through the charge states that are evenly divisible and exclude this cluster if it shares core peaks with a previously called cluster
                foreach (var chargeFactorCluster in (finalResult.Where(c => (c.Z > aCluster.Z) && (c.Z % aCluster.Z == 0))))
                {
                    // Look for *significant* overlap between the core peaks of these two ion clusters

                    int skipMin = (int)(Math.Min(chargeFactorCluster.Peaks.Count(p => p.IsCorePeak), aClusterCorePeakCount) * 0.35);  // core peaks overlap can be no more than 30%
                    int skipPoints = 0;

                    foreach (var aPeak in chargeFactorCluster.Peaks.Where(y => y.IsCorePeak))
                    {
                        if (aCluster.Peaks.Where(p => p.IsCorePeak).Any(x => x.MZ == aPeak.MZ))
                        {
                            skipPoints++;
                        }
                    }

                    if (skipPoints > skipMin)
                    {
                        skip = true;
                        break;
                    }
                }

                if (!skip) finalResult.Add(aCluster);

            }

        }

        var unAssignedMzList = new SortedList<double, float>(mzList);
        var avgCoreIntensity = 0d;
        var coreCount = 0;

        foreach (var aCluster in finalResult)
        {
            foreach (var aPeak in aCluster.Peaks)
            {
                if (aPeak.IsCorePeak)
                {
                    coreCount++;
                    avgCoreIntensity += aPeak.Intensity;
                }

                if (mzList.ContainsKey(aPeak.MZ)) unAssignedMzList.Remove(aPeak.MZ);

            }
        }

        avgCoreIntensity = avgCoreIntensity / coreCount;

        var t = unAssignedMzList.Where(p => p.Value > avgCoreIntensity);
        Cluster currentCluster = null;

        foreach (var aCharge in AssumedCharges)
        {
            foreach (var peak in t)
            {
                currentCluster = new Cluster() { Peaks = new ListOfClusterPeak(), Z = aCharge, _monoMZ = peak.Key, Score = 301 };
                currentCluster.Peaks.Add(new ClusterPeak() { MZ = peak.Key, Z = aCharge, Index = 0, IsCorePeak = true, Intensity = peak.Value });
                finalResult.Add(currentCluster);
            }
        }

        return finalResult;

        //TODO: Filter out charge states that are a factor of another...
    }

    public static int[] AssumedCharges = new int[] { 1 };


    /// <summary>
    /// Remove peaks that are not part of a group and remove peaks with groups that are missing major peaks
    /// </summary>
    /// <param name="result"></param>
    /// <param name="zSpace"></param>
    /// <returns></returns>
    private List<Cluster> FindIsotopeGroups(List<double> input, int charge, double massTolerance, bool applyNoiseFilter = true)
    {
        // declare vars
        var clusters = new List<Cluster>();
        double LooseTolerance = 0.01;
        double dCharge = (double)charge;
        var zSpace = SignalProcessor.ProtonMass / dCharge;  // changed to accomodate mass defect

        foreach (var aPeak in input)
        {
            //TODO: use cross-correlation here to find peaks of the distribution?

            if (clusters.Any(c => c.Peaks.Any(p => p.MZ == aPeak && c.Peaks.Any(x => x.IsMonoisotopic) && p.Index >= c.Peaks.Where(x => x.IsMonoisotopic).First().Index)))
                continue;

            LooseTolerance = ((0.06 * massTolerance) / 1e6d) * (aPeak * dCharge);

            LooseTolerance = Math.Max(LooseTolerance, zSpace / 5d);

            var currentCluster = new Cluster() { Z = charge };
            
            var distribution = new System.Collections.Generic.Dictionary<int, ClusterPeak>();

            double peakMass = aPeak * dCharge * SignalProcessor.ProtonMass;

            // Derived from polynomial fit to a in-silico computed plot of averagine for first/last peak from apex of averagine
            int peakSearchStart = (int)(-0.000000004 * (peakMass * peakMass) + 0.0007 * peakMass - 0.8035);
            peakSearchStart = Math.Max(4, (peakSearchStart + 2)) * -1;  // old and good
            int peakSearchEnd = (int)(-0.00000002 * (peakMass * peakMass) + 0.001 * peakMass + 3.6539);
            peakSearchEnd = Math.Max(4, peakSearchEnd + 2);  //  I add 2 in case we miss the apex by up to 2 peaks.  

            int peakSearchCount = peakSearchEnd - peakSearchStart;

            int firstRunPeak = 0;
            int lastRunPeak = 0;
            int runLength = 0;
            bool inAnchorRun = false;
            bool foundAnchorRun = false;

            int currentRunStart = 0;

            for (int i = 0; i < peakSearchCount; i++)
            {
                var mzTarget = aPeak + ((i + peakSearchStart) * zSpace);
                var candidate = mzList.FindClosest(mzTarget);

                var newPeak = new ClusterPeak() { Parent = (currentCluster), MZ = candidate, Intensity = mzList[candidate], Index = i, Delta = (float)Math.Abs(mzTarget - candidate) };

                if (newPeak.Delta > LooseTolerance)
                {
                    // Could not find a peak with a close enough match, so we are throwing out the match and using a 0-intensity dummy peak
                    if (newPeak.Delta > (LooseTolerance))
                        newPeak = new ClusterPeak() { Parent = (currentCluster), MZ = mzTarget, Intensity = 0, Index = i, Delta = float.MaxValue };

                    if (!foundAnchorRun)
                    {
                        currentRunStart = newPeak.Index;
                        if (inAnchorRun) foundAnchorRun = true;
                    }
                }
                else if (runLength < (i - currentRunStart) + 1)
                {
                     var minIntensityPercent = 0.1509 * Math.Log(newPeak.Mass) - 0.6479;

                    minIntensityPercent = minIntensityPercent * 0.5;  //  Add some tolerance factor (it has to be within 30% of expected value for neighboring intensity delta)?  This should be forgiving because ion suppression may cause real variability in neighboring peaks?

                    if (i > 0 && aPeak != mzTarget && (Math.Min(newPeak.Intensity, distribution[i - 1].Intensity) / Math.Max(newPeak.Intensity, distribution[i - 1].Intensity) < minIntensityPercent))
                    {
                        currentRunStart = i;
                        if (!foundAnchorRun && inAnchorRun) foundAnchorRun = true;
                    }
                    else
                    {
                        runLength = (i - currentRunStart) + 1;
                        firstRunPeak = currentRunStart;
                        lastRunPeak = i;
                    }
                }

                // This run contains the anchor peak! 
                if (aPeak == mzTarget)
                    inAnchorRun = true;

                 distribution.Add(i, newPeak);
            }

            // Ion spacing scoring
            int score = 0;
            if (runLength <= 1) continue;  // abort if we can't find a decent length run of peaks

            //MDK  changes 
            // ion run too small
            if (charge < 3 && runLength < 2) continue;
            if (charge > 5 && runLength < 4) continue;
            if (charge > 6 && runLength < 6) continue;
            if (charge > 12 && runLength < 9) continue;

            // ion run to long
            if (charge < 3 && runLength > 5) continue;
            if (charge < 5 && runLength > 14) continue;

            score += 10 * runLength;

            double VarianceInRun = 0;
            double AverageDeviation = 0;
            double localNoise = 0;

            float secondLargestIntensity = 0;

            try
            {
                secondLargestIntensity = distribution.Where(p => p.Key >= firstRunPeak && p.Key <= lastRunPeak).OrderByDescending(k => k.Value.Intensity).Skip(1).Take(1).First().Value.Intensity;
            }
            catch
            {
                // This should be very rare.  There is no second most intense peak...
                continue;
            }

            currentCluster.Peaks = new ListOfClusterPeak(distribution.Values);
            currentCluster.SetMonoAndCorePeaks();

            // Find the Core peak group that contains the initial peak (aPeak)

            var mainPeak = currentCluster.Peaks.Where(p => p.MZ == aPeak).First();
            if (currentCluster.Peaks.Count(p => p.IsCorePeak) == 1) currentCluster.Peaks[mainPeak.Index + 1].IsCorePeak = true;  // This will make sure that we have the minimum 2 core peaks.  However, we my not have a correct intensity pattern.  

            int coreStartIndex = mainPeak.Index;
            int coreEndIndex = mainPeak.Index;

            foreach (var currentPeak in currentCluster.Peaks.Where(p => p.Index < mainPeak.Index).OrderByDescending(x => x.Index))
            {
                if (!currentPeak.IsCorePeak)
                {
                    coreStartIndex = currentPeak.Index + 1;
                    break;
                }
            }

            foreach (var currentPeak in currentCluster.Peaks.Where(p => p.Index > mainPeak.Index).OrderBy(x => x.Index))
            {
                if (!currentPeak.IsCorePeak)
                {
                    coreEndIndex = currentPeak.Index - 1;
                    break;
                }
            }

            ClusterPeak lastPeak = null;

            // Start with the coreStartIndex, then calculate the stats for just the true core!  

            foreach (var currentPeak in currentCluster.Peaks.Where(p => p.IsCorePeak && p.Index >= coreStartIndex && p.Index <= coreEndIndex))
            {
                
                if (lastPeak == null)
                {
                    lastPeak = currentPeak;
                    continue;
                }
                
                AverageDeviation += currentPeak.MZ - lastPeak.MZ;
                VarianceInRun += ((currentPeak.MZ - lastPeak.MZ) - (SignalProcessor.ProtonMass / (double)charge)) * ((currentPeak.MZ - lastPeak.MZ) - (SignalProcessor.ProtonMass / (double)charge));

                for (int j = mzList.FindClosestIndex(lastPeak.MZ) + 1; j < mzList.FindClosestIndex(currentPeak.MZ); j++)
                {
                    localNoise += (mzList.Values[j] > ((Math.Min(lastPeak.Intensity, currentPeak.Intensity)) * 0.65)) ? Math.Min(1, mzList.Values[j] / secondLargestIntensity) : 0;
                }
                if (verbose) Debug.Print("Current Peak:  " + currentPeak.MZ);
                if (verbose) Debug.Print("Local Noise:  " + localNoise);
                lastPeak = currentPeak;
            }

            if (verbose) Debug.Print("Evaluating Distribution of length = " + runLength.ToString() + ", charge = " + charge.ToString());

            
            if (applyNoiseFilter && (localNoise / ((coreEndIndex - coreStartIndex) * zSpace)) > 4)
            {
                if (verbose) Debug.Print("Reject - Local Noise " + (localNoise / (coreEndIndex - coreStartIndex)).ToString("0.0000"));
                continue;
            }
            // DEBUG

            // WRONG, this is not the core count!!!!

            VarianceInRun /= coreEndIndex - coreStartIndex;
            AverageDeviation /= coreEndIndex - coreStartIndex;
            var StandardDev = Math.Sqrt(VarianceInRun);

            // Detect Charge from Run and compare...
            var detectedCharge = SignalProcessor.ProtonMass / AverageDeviation;  // Calculate Charge state from the run of peaks 

            if (Math.Abs(detectedCharge - (double)charge) > 0.25)
            {
                if (verbose) Debug.Print("Reject - Charge State Validate Failure.  " + detectedCharge.ToString() + " detected, " + charge.ToString() + " expected");
                continue;  // Charge state Validation Failure
            }

            score += (int)(100d / ((Math.Abs(detectedCharge - (double)charge) * 8d) + 0.01));
            score += (int)(0.03 / (StandardDev + 0.000001));

            currentCluster.Score = score;
            if (verbose) Debug.Print("Current Cluster Score: " + score.ToString() );

            // Calling FindMono here because the next "if" statement would call it automatically by checking the "IsMonoisotopic property, but we want to prepopulate, so we don't try to populate later
            //currentCluster.FindMono();
            //currentCluster.FindMonoMz();

            //currentCluster.FindMonoMzWithMin();

            //currentCluster.FindCorePeaks();

            //currentCluster.SetMonoAndCorePeaks();

            if (!currentCluster.Peaks.Where(p => p.IsCorePeak).Any())
            {
                if (verbose) Debug.Print("Rejected: No corePeaks detected!");
                continue;
            }

            score += currentCluster.Peaks.Where(p => p.IsCorePeak).Count() * 30;


            //MDK increase thershold to 1000 for score on Agilent Data
            if (score < 1000) continue;

            if (score < 300 || (currentCluster.Peaks[0].Mass < 3000 && !currentCluster.Peaks.Any(p => p.IsMonoisotopic && p.IsCorePeak)))
            {
                if (verbose && score > 300)
                    Debug.Print("Rejected because all charge 3 and under ions must have a monoisotopic peak and it must be a core peak");
                else if (verbose)
                    Debug.Print("Rejected for low score = " + score.ToString());

                continue;
            }

            // Estimate local noise level by intensity of peaks between isotope peaks...

            if (verbose) Debug.Print("Score = " + score.ToString() + "; A good one at charge " + charge.ToString() + "?  " + runLength.ToString() + " continuous peaks at m/z => " + distribution[firstRunPeak].MZ.ToString("0.0000") + " to " + distribution[lastRunPeak].MZ.ToString("0.0000"));

            clusters.Add(currentCluster);
                                
          }
        
        return clusters;
    }
  

    public List<Cluster> ParallelDetectChargeStates(int MinCharge, int MaxCharge, double massTolerance)
    {
        // m/z, charge state, cluster ID
        //var clusters = new MultiKeyDictionary<double, int, int>();

        // Can be Parallel!!!
        var finalResult = new List<Cluster>();

        //foreach (var chargeState in Enumerable.Range(MinCharge, (MaxCharge - MinCharge) + 1).Reverse())

        Parallel.ForEach((Enumerable.Range(MinCharge, (MaxCharge - MinCharge) + 1).Reverse()), chargeState =>
        {
            var result = DetectCharge(chargeState, massTolerance);


            //TODO: filter out peaks that don't have others from the isotope distribution
            //TODO: filter out distributions that have gaps (or missing peaks) -- open up tolerance a bit to be sure though...

            //double zSpace = 1.007947d / (double)chargeState;

            // Remove peaks with insufficient supporting peak data (i.e. at least 3 peaks where there should be 6)
            //result.RemoveWhere(t => result.Where(x => ((x < (t + (3.2 * zSpace))) && (x > (t - (3.2 * zSpace))))).Count() < 4);

            //TODO: take the groupings of peaks and fill in the ones missing from the distribution.  
            //result

            var result2 = FindIsotopeGroups(result, chargeState, massTolerance);

            foreach (var aCluster in result2)
                finalResult.Add(aCluster);

            //                finalResult.Add(chargeState, result);

            //foreach (var aPeak in result)
            //    Debug.WriteLine("Found a peak (Z=" + chargeState + "): " + aPeak);

            //Debug.Print("Charge States Assigned = " + result.Count);



        });

        return finalResult;

        //TODO: Filter out charge states that are a factor of another...
    }

    public List<double> DetectCharge(int state, double massTolerance, double massRangeMin = 0, double massRangeMax = double.MaxValue)
    {
        //Debug.WriteLine("Detecting Charge for " + state.ToString());

        //double target
        var ChargeDetected = new List<double>();

        //double chargeSpace = 1.007947d / (double)state;
        double chargeSpace = SignalProcessor.ProtonMass / (double)state;  // changed to accomodate mass defect

        double hdt = 0.05 + (0.4f * Math.Min(((float)state / 17f), 1f)); // Height difference tolerance  -- multiplier that varies based on charge state.  
        // At Higher Charge, we are less tolerant to adjacent peaks having an intensity difference.  

        var toleranceBase = massTolerance / 1e6d;  // accessing the settings is an expensive operation, so let's do it outside the loop!

        double tolerance = 0;

        // convert mass values to a charge specific m/z
        var mzMin = massRangeMin.ToMZ(state);
        var mzMax = massRangeMax == double.MaxValue ? double.MaxValue : massRangeMax.ToMZ(state);  // the MaxValue check is to prevent an overflow when calculating the Mass

        // Reverse for loop (forr + tab)
        // Take all peaks
        for (int i = intensityList.Count - 1; i >= 0; --i)
        // Take top 80% peaks
        //for (int i = intensityList.Count - 1; i >= (int)(intensityList.Count * 0.2); --i)  // take only top 80% most intense peaks -- arbitrary
        
        {
            double aroundMZ = intensityList[intensityList.Keys[i]];

            if (aroundMZ < mzMin || aroundMZ > mzMax) continue;  // skip peaks that are outside the target mass range (for better performance when used)

            double mz = mzList.FindClosest(aroundMZ + chargeSpace);

            tolerance = toleranceBase * mz;

            //var tolerance_old = Tolerance * (state > 2 ? 1 : 1.7);  // more loose tolerance for the lower charge ions

            if ((Math.Abs((aroundMZ + chargeSpace) - mz) < tolerance) && (Math.Min(mzList[mz], mzList[aroundMZ]) > Math.Max(mzList[mz], mzList[aroundMZ]) * hdt))
            {
                if (state > 2)
                {
                    double mz2 = mzList.FindClosest(aroundMZ - chargeSpace);
                    if ((Math.Abs((aroundMZ - chargeSpace) - mz2) < tolerance) && (Math.Min(mzList[mz2], mzList[aroundMZ]) > Math.Max(mzList[mz2], mzList[aroundMZ]) * hdt)) ChargeDetected.Add(aroundMZ);
                }
                else
                {
                    ChargeDetected.Add(aroundMZ);
                }
            }
         }
        //Debug.WriteLine("*** " + ChargeDetected.Count.ToString() + " ion spaces detected for charge " + state.ToString());
        return ChargeDetected;
    }
}


public static class MassSpec
{
    public static double MonoToApexMass(double mass)
    {
        var theoIon = GenerateAveragine(mass);

        var pattern = MassSpectrometry.IsotopeCalc.CalcIsotopePeaks(theoIon, 1f);

        string resultStr = string.Empty;

        foreach (var peak in pattern)
            resultStr += peak.Key.ToString() + ", " + peak.Value.ToString() + "\n";

        return pattern.MaxBy(p => p.Value).Key - pattern.Keys[0];

    }

    // Generates approximate chemical formula for therotical isotope pattern display
    internal static Ion GenerateAveragine(double mw)
    {
        //MDK Determine what isotope table is loaded
        AveragineCacheSettings settings = AveragineCacheSettings.Instance;
        var newCacheName = settings.SelectedCacheFile;
        var formula = "";

        //For Proteins and Peptides
        if (newCacheName.Contains("Protein"))
            {
            float AveragineMW = 111.1254f;
            float AverageC = 4.9384f;
            float AverageH = 7.7583f;
            float AverageN = 1.3577f;
            float AverageO = 1.4773f;
            float AverageS = 0.0417f;

            var roundNumAveragine = (int)Math.Round(mw / AveragineMW, 0);

            // Example: C(644) H(1012) N(177) O(193) S(5)
            var proteinFormula = "C(" + Math.Round(AverageC * roundNumAveragine, 0)
                      + ") H(" + Math.Round(AverageH * roundNumAveragine, 0)
                      + ") N(" + Math.Round(AverageN * roundNumAveragine, 0)
                      + ") O(" + Math.Round(AverageO * roundNumAveragine, 0)
                      + ") S(" + Math.Round(AverageS * roundNumAveragine, 0) + ")";

            formula = proteinFormula;
            }

        // for oligos
        if (newCacheName.Contains("Oligo"))
            {
            float AveragineMW = 305.8335f;
            float AverageC = 9.75f;
            float AverageH = 12.30f;
            float AverageN = 3.75f;
            float AverageO = 5.90f;
            float AverageP = 0.95f;

            var roundNumAveragine = (int)Math.Round(mw / AveragineMW, 0);

            var oligoFormula = "C(" + Math.Round(AverageC * roundNumAveragine, 0)
                  + ") H(" + Math.Round(AverageH * roundNumAveragine, 0)
                  + ") N(" + Math.Round(AverageN * roundNumAveragine, 0)
                  + ") O(" + Math.Round(AverageO * roundNumAveragine, 0)
                  + ") P(" + Math.Round(AverageP * roundNumAveragine, 0) + ")";

            formula = oligoFormula ;
            }

        return new Ion(formula, -1);
    }
}

public static partial class Extensions
{
    /// <summary>
    /// Find the closest X value to the target passed in using a binary search.  
    /// </summary>
    /// <param name="target">Target X Value</param>
    /// <returns>Closest X Value or double.NaN if specified tolerance is not met</returns>
    public static double FindClosest(this SortedList<double, float> list, double target, double tolerance = -1)
    {
        // http://www.brpreiss.com/books/opus6/html/page452.html


        // if (list.ContainsKey(target)) return target;

        int left = 0, middle = 0, right = list.Count - 1;

        while (left <= right)
        {
            middle = (left + right) / 2;
            if (target > list.Keys[middle])
                left = middle + 1;
            else if (target < list.Keys[middle])
                right = middle - 1;
            else
                return (tolerance < 0 || (Math.Abs(list.Keys[middle] - target) < tolerance)) ? list.Keys[middle] : double.NaN;
        }

        // Validate bounds!
        left = Math.Min(left, (list.Keys.Count - 1));
        right = Math.Max(right, 0);

        if (Math.Abs(target - list.Keys[left]) > Math.Abs(target - list.Keys[right]))
            return (tolerance < 0 || (Math.Abs(list.Keys[right] - target) < tolerance)) ? list.Keys[right] : double.NaN;
        else
            return (tolerance < 0 || (Math.Abs(list.Keys[left] - target) < tolerance)) ? list.Keys[left] : double.NaN;
    }

    /// <summary>
    /// Find the closest X value to the target passed in using a binary search.  
    /// </summary>
    /// <param name="target">Target X Value</param>
    /// <returns>Closest X Value or double.NaN if specified tolerance is not met</returns>
    public static int FindClosestIndex(this SortedList<double, float> list, double target, double tolerance = -1)
    {
        // http://www.brpreiss.com/books/opus6/html/page452.html


        // if (list.ContainsKey(target)) return target;

        int left = 0, middle = 0, right = list.Count - 1;

        while (left <= right)
        {
            middle = (left + right) / 2;
            if (target > list.Keys[middle])
                left = middle + 1;
            else if (target < list.Keys[middle])
                right = middle - 1;
            else
                return (tolerance < 0 || (Math.Abs(list.Keys[middle] - target) < tolerance)) ? middle : -1;
        }

        // Validate bounds!
        left = Math.Min(left, (list.Keys.Count - 1));
        right = Math.Max(right, 0);

        if (Math.Abs(target - list.Keys[left]) > Math.Abs(target - list.Keys[right]))
            return (tolerance < 0 || (Math.Abs(list.Keys[right] - target) < tolerance)) ? right : -1;
        else
            return (tolerance < 0 || (Math.Abs(list.Keys[left] - target) < tolerance)) ? left : -1;
    }

    private static double GetSigma()
    {
        return GetSigma(30000);
    }

    private static double GetSigma(int resolution)
    {
        double m = 1000;  // value suggested by Rick and Mike

        double fwhm = m / resolution;

        // Based on: http://mathworld.wolfram.com/GaussianFunction.html
        return fwhm / 2.35482;  // (2 * Math.Sqrt(2 * Math.Log(2)));
    }

    private static double[] GenerateGaussian()
    {
        return GenerateGaussian(30000, 0.05f);
    }

    private static double[] GenerateGaussian(int resolution, float binSize)
    {
        var points = new List<double>();

        var s = GetSigma(resolution);

        //System.Diagnostics.Debug.WriteLine("Gaussian internal");

        for (int i = -10; i <= 10; i++)
        {
            var x = i * binSize;

            var g = (Math.Exp((x * x) / (-2 * s * s)) / (s * Math.Sqrt(2 * Math.PI)));

            points.Add(g);

            //System.Diagnostics.Debug.WriteLine(x.ToString() + "\t " + g.ToString());

        }

        //System.Diagnostics.Debug.WriteLine(" End Gaussian internal");

        return points.ToArray();
    }

    private static int IndexOfMax(double[] result)
    {
        int maxIndex = 0;

        for (int i = 0; i < result.Length; i++)
            if (result[i] > result[maxIndex]) maxIndex = i;

        return maxIndex;
    }

    public class point
    {
        //    [OperationContract]
        public double X { get; set; }
        //    [OperationContract]
        public float Y { get; set; }

        public override string ToString()
        {
            return X.ToString() + ", " + Y.ToString();
        }
    }


}


