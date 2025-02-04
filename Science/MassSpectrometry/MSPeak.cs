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

// Based on Project Morpheus
// https://github.com/cwenger/Morpheus/blob/master/Morpheus/mzML/TandemMassSpectra.mzML.cs

namespace MassSpectrometry
{
    public class MSPeak
    {
        public double MZ { get; private set; }

        public double Intensity { get; private set; }

        public int Charge { get; private set; }

        public double Mass { get; private set; }

        public MSPeak(double mz, double intensity, int charge, int polarity)
        {
            MZ = mz;
            Intensity = intensity;
            Charge = charge;
            CalculateMass(charge, polarity);
        }

        private void CalculateMass(int charge, int polarity)
        {
            if(charge == 0)
            {
                charge = polarity;
            }
            Mass = MassFromMZ(MZ, charge);
        }

        public static double MassFromMZ(double mz, int charge)
        {
            return charge == 0 ? mz : mz * Math.Abs(charge) - charge * Constants.PROTON_MASS;
        }

        public static double MZFromMass(double mass, int charge)
        {
            if(charge == 0)
            {
                throw new ArgumentOutOfRangeException("Charge cannot be zero.");
            }
            return (mass + charge * Constants.PROTON_MASS) / Math.Abs(charge);
        }

        public static int AscendingMZComparison(MSPeak left, MSPeak right)
        {
            return left.MZ.CompareTo(right.MZ);
        }

        public static int DescendingIntensityComparison(MSPeak left, MSPeak right)
        {
            return -(left.Intensity.CompareTo(right.Intensity));
        }

        public static int AscendingMassComparison(MSPeak left, MSPeak right)
        {
            return left.Mass.CompareTo(right.Mass);
        }
    }
}
