package com.example.sun;

import java.text.DecimalFormat;
import java.util.Scanner;

public class MoonCalculator {

    //moon
    public static double eccentricAnomaly = 0;

    public static final double ECCENTRICITYMOON = 0.0549;

    public static final double ECLIPTICLONGMOONEPOCH = 218.316433;//ecliptic longitude of the moon at the epoch (degrees)
    public static double meanEclipticLongMoon = 0;//mean ecliptic longitude of the moon (degrees)
    public static double meanEclipticLongMoonAscending = 0;//mean ecliptic longitude of the moon ascending node(degrees)
    public static double meanAnomalyMoon = 0;//mean anomaly of the moon, uncorrected
    public static double annualEquationCorrection = 0;//annual equation correction
    public static double evectionCorrection = 0;//evection correction
    public static double meanAnomalyCorrection = 0;//mean anomaly correction
    public static final double ANGSPEEDMOON = 13.176339686;//angular speed of moon along its mean orbit (degrees per day)
    public static final double ECLIPTICLONGMOONASCENDING = 125.044522;//ecliptic longitude of the moon ascending node at the epoch (degrees)
    public static final double ECLIPTICLONGMOONASCENDINGCONST = 0.0529539;//constant/coefficient to calculate ecliptic longitude of the moon ascending node at the epoch
    public static final double MEANANOMALYCONST = 0.111404;//constant/coefficient to calculate mean anomaly of the moon, uncorrected
    public static final double ECLIPTICLONGMOONPERIGEE = 83.353451;//ecliptic longitude of the moon at perigee at the epoch (degrees)
    public static final double ANNUALEQUATIONCONST = 0.1858;//constant/coefficient to calculate annual equation
    public static final double EVECTIONCONST = 1.2739;//constant/coefficient to calculate evection correction
    public static final double MEANANOMCORRECTIONCONST = 0.37;//constant/coefficient to calculate mean anomaly correction

    //calculate moon's mean ecliptic longitude
    public double meanEclipticLong(double days) {
        double lambda = (ANGSPEEDMOON * days) + ECLIPTICLONGMOONEPOCH;
        meanEclipticLongMoon = lambda%360;
        return meanEclipticLongMoon;
    }

    //calculate moon's mean ecliptic longitude for the ascending node
    public double meanEclipticLongAscendingNode(double days) {
        double omega = ECLIPTICLONGMOONASCENDING - (ECLIPTICLONGMOONASCENDINGCONST * days);
        meanEclipticLongMoonAscending = omega%360;
        if (meanEclipticLongMoonAscending < 0) {
            meanEclipticLongMoonAscending = meanEclipticLongMoonAscending + 360;
        }
        return meanEclipticLongMoonAscending;
    }

    //calculate moon's mean anomaly
    public double meanAnomaly(double days) {
        double mAnom = meanEclipticLongMoon - (MEANANOMALYCONST * days) - ECLIPTICLONGMOONPERIGEE;
        meanAnomalyMoon = mAnom%360;
        if (meanAnomalyMoon < 0) {
            meanAnomalyMoon = meanAnomalyMoon + 360;
        }
        return meanAnomalyMoon;
    }

    //calculate annual equation correction
    public double annualEquation(double meanAnomalySun) {
        annualEquationCorrection = ANNUALEQUATIONCONST * (Math.sin(degreesToRadians(meanAnomalySun)));
        return annualEquationCorrection;
    }

    //calculate evection correction
    public double evection(double eclipticLongitudeSun) {
        evectionCorrection = EVECTIONCONST * (Math.sin(degreesToRadians((2*(meanEclipticLongMoon-eclipticLongitudeSun)-meanAnomalyMoon))));
        return evectionCorrection;
    }

    //calculate mean anomaly correction
    public double meanAnomCorrection(double meanAnomalySun) {
        meanAnomalyCorrection = meanAnomalyMoon + evectionCorrection - (MEANANOMCORRECTIONCONST*(Math.sin(degreesToRadians(meanAnomalySun))));
        return meanAnomalyCorrection;
    }

    //convert degrees to radians
    public static double degreesToRadians(double degrees)
    {
        return degrees * 0.0174532925199;
    }

    //convert radians to degrees
    public static double radiansToDegrees(double radians)
    {
        return radians / 0.0174532925199;
    }
}
