package com.example.sun;

import java.text.DecimalFormat;
import java.util.Scanner;

public class SunCalcOld {

    public static Scanner keyboard = new Scanner(System.in);//input stream

    private static DecimalFormat df2 = new DecimalFormat("00.00");
    private static DecimalFormat dfw = new DecimalFormat("##");
    private static DecimalFormat dfm = new DecimalFormat("00");

    public static double year = 0;
    public static double month = 0;
    public static double day = 0;
    public static double timeDay = 0;
    public static double julianDayRaw = 0;//Julian Day number not accounting for time of day. (Ends .5 - exact time at midnight)
    public static double julianDay = 0;//Julian Day number exact
    public static double meanAnomDeg = 0;//the mean anomaly in degrees
    public static double trueAnomaly = 0;//true anomaly
    public static double eclipticLongSun = 0;//ecliptic longitude of the sun
    public static double rASun = 0;//right ascension of sun
    public static double decSun = 0;//declination of sun
    public static double latitude = 0;//earth latitude, north positive   54.6 (nl 52)
    public static double longitude = 0;//earth longitude, east positive   -5.9 (nl 5)
    public static double jTransit = 0; //j (time of transit, measured in terms of Julian date)
    public static double jRise = 0;//sunrise time (Julian date)
    public static double jSet = 0;//sunset time (Julian date)
    public static double hourAngle = 0;

    public static final double J2K = 2451545;//J2000
    public static final double ELOP = 102.9373;//ecliptic longitude of perihelion
    public static final double OOTE = 23.4393;//obliquity of the ecliptic
    public static final double JZERO = 0.0009;
    public static final double JONE = 0.0053;
    public static final double JTWO = -0.0068;
    public static final double M0 = 357.5291;//mean anomaly degrees
    public static final double M1 = 0.98560028;//mean anomaly degrees per day
    public static final double THETAZERO = 280.1470;//sidereal constant for earth
    public static final double THETAONE = 360.9856235;//sidereal constant for earth

    //main method
    public static void main(String[] args) {
        getDate();
        getTime();
        getPlace();
        dateTimeToJulianDay();
        System.out.println("LATITUDE/LONGITUDE ------------------------------------------------------------ " + latitude + ", " + (-longitude));
        meanAnomaly();
        equationCentre();
        eclipticalCoordinates();
        equatorialCoordinates();
        observer();
        solarTransit();//repeat stages for calculated value of jTransit to increase accuracy ~ meanAnomalyFromJulianDayTime(jTransit);
        solarTransitTime();
        riseAndSet();
        riseAndSetTime();
    }

    //get calender date from user
    public static void getDate() {
        System.out.println("Enter year: ");
        year = keyboard.nextDouble();
        System.out.println("Enter month: ");
        month = keyboard.nextDouble();
        System.out.println("Enter day: ");
        day = keyboard.nextDouble();
    }

    //get time of day from user
    public static void getTime() {
        System.out.println("Enter time of day:(24hh.mm) ");
        timeDay = keyboard.nextDouble();
    }

    //get latitude and longitude from user
    public static void getPlace() {
        System.out.println("Enter latitude (north positive): ");
        latitude = keyboard.nextDouble();
        System.out.println("Enter longitude (east positive): ");
        longitude = -(keyboard.nextDouble());
    }

    //convert (gregorian) calender date to julian day number - exact method including time of day
    //exact whole number date gives 12.00 hours (noon), date given as e.g. 0.5 corresponds to 00.00 hours (midnight, UTC)
    public static double dateTimeToJulianDay() {
        double hour = Math.floor(timeDay);
        double minutes = timeDay - hour;
        double minDecFraction = (minutes/60)*100;
        double decTime = (hour + minDecFraction)/24;

        //keep year and month as entered for display
        double yearDisplay = year;
        double monthDisplay = month;

        if ((month == 1) || (month == 2)) {
            year = (year - 1);
            month = (month + 12);
        }

        double a = Math.floor(year/100);
        double b = Math.floor(a/4);
        double c = 2 - a + b;
        double e = Math.floor(365.25*(year+4716));
        double f = Math.floor(30.6001*(month+1));

        julianDayRaw = c + day + e + f - 1524.5;

        System.out.println("DATE/TIME --------------------------------------------------------------------- " + dfw.format(yearDisplay) + "/"
                + dfw.format(monthDisplay) + "/" + dfw.format(day) + " " + df2.format(timeDay) + " (UTC)");
        System.out.println("JULIAN DAY (DATE ONLY) -------------------------------------------------------- " + julianDayRaw);

        julianDay = julianDayRaw + decTime;
        System.out.println("EXACT JULIAN DAY NUMBER (INCLUDING TIME) -------------------------------------- " + julianDay);
        return julianDay;
    }

    //calculate the mean anomaly in degrees
    public static double meanAnomaly() {
        double theta= M0 + M1*(julianDay - J2K);
        meanAnomDeg = theta%360;

        System.out.println("MEAN ANOMALY, EARTH (DEGREES) ------------------------------------------------- " + meanAnomDeg);
        return meanAnomDeg;
    }

    //calculate the mean anomaly in degrees
    public static double meanAnomalyFromJulianDayTime(double julianDayTime) {
        double theta= M0 + M1*(julianDayTime - J2K);
        meanAnomDeg = theta%360;

        System.out.println("TRANSIT MEAN ANOMALY, EARTH (DEGREES) ----------------------------------------- " + meanAnomDeg);
        return meanAnomDeg;
    }

    //calculate the equation of centre - trig functions take arguments in radians
    public static double equationCentre() {
        final double C1 = 1.9148;
        final double C2 = 0.0200;
        final double C3 = 0.0003;

        //convert to radians
        double meanAnomRad = degreesToRadians(meanAnomDeg);
        double equationOfCentre = (C1 * Math.sin(meanAnomRad)) + (C2 * Math.sin(2 * meanAnomRad)) + (C3 * Math.sin(3 * meanAnomRad));
        trueAnomaly = meanAnomDeg + equationOfCentre;

        System.out.println("EQUATION OF CENTRE ------------------------------------------------------------ " + equationOfCentre);
        System.out.println("TRUE ANOMALY ------------------------------------------------------------------ " + trueAnomaly);
        return trueAnomaly;
    }

    //calculates the ecliptic longitude of the sun
    public static double eclipticalCoordinates() {
        double lambda = trueAnomaly + ELOP + 180;
        eclipticLongSun = lambda%360;

        System.out.println("ECLIPTIC LONGITUDE OF SUN ----------------------------------------------------- " + eclipticLongSun);
        return eclipticLongSun;
    }

    //calculates the right ascension and declination of the sun
    public static double equatorialCoordinates() {
        //constant approximations for right ascension and declination (earth)
        //not needed using alternative method only

        //convert to radians
        double eclipticLongSunRad = degreesToRadians(eclipticLongSun);
        double sinLambdaSun = Math.sin(eclipticLongSunRad); //sine of the ecliptic longitude
        double cosLambdaSun = Math.cos(eclipticLongSunRad); //cosine of the ecliptic longitude
        //calculate right ascension using alternative method
        //??sun = arctan(sin ??sun cos ??, cos ??sun). Returns angle in radians.
        rASun = Math.atan2(sinLambdaSun*Math.cos(degreesToRadians(OOTE)), cosLambdaSun);
        rASun = radiansToDegrees(rASun);
        if (rASun < 0) {
            rASun = 360 + rASun; //convert to range 0-360 if negative. What about greater than or equal to 360?
        }
        //calculate declination using alternative method
        //??sun = arcsin(sin ??sun sin ??). Returns angle in radians.
        decSun = Math.asin(sinLambdaSun*Math.sin(degreesToRadians(OOTE)));
        decSun = radiansToDegrees(decSun);

        System.out.println("SUN'S RIGHT ASCENSION (ALTERNATIVE METHOD) --------- " + rASun);
        System.out.println("SUN'S DECLINATION (ALTERNATIVE METHOD) ------------- " + decSun);
        return rASun;
        //are alternative method values more accurate? Practically the same. Assume full equations are slightly better.
    }

    //calculates sun position for earth observer
    public static double observer() {
        double thetaSidereal = THETAZERO + THETAONE * (julianDay - J2K) - longitude; //raw angle, remove multiples of 360
        double sidereal = thetaSidereal%360; //sidereal time
        hourAngle = sidereal - rASun;

        System.out.println("SIDEREAL TIME (DEGREES) ------------------------------------------------------- " + sidereal);
        System.out.println("HOUR ANGLE H ------------------------------------------------------------------ " + hourAngle);

        //convert to radians
        double hourAngleRad = degreesToRadians(hourAngle);
        double latitudeRad = degreesToRadians(latitude);
        double decSunRad = degreesToRadians(decSun);
        double azimuth = Math.atan2(Math.sin(hourAngleRad), Math.cos(hourAngleRad)*Math.sin(latitudeRad) - Math.tan(decSunRad)*Math.cos(latitudeRad));
        azimuth = radiansToDegrees(azimuth) + 180; //convert from radians to degrees, add 180 degrees for conventional azimuth (north=0)

        double altitude = Math.asin(Math.sin(latitudeRad)*Math.sin(decSunRad) + Math.cos(latitudeRad)*Math.cos(decSunRad)* Math.cos(hourAngleRad));
        altitude = radiansToDegrees(altitude); //convert from radians to degrees
        System.out.println("AZIMUTH ----------------------------------------------------------------------- " + azimuth);
        System.out.println("ALTITUDE ---------------------------------------------------------------------- " + altitude);
        //azimuth 0 = north
        return sidereal;
    }

    public static double hourAngleFromJulianDayTime(double julianDayTime) {
        double thetaSidereal = THETAZERO + THETAONE * (julianDayTime - J2K) - longitude; //raw angle, remove multiples of 360
        double sidereal = thetaSidereal%360; //sidereal time
        hourAngle = sidereal - rASun;

        return hourAngle;
    }

    //calculate transit time (julian day)
    public static double solarTransit() {
        double n;
        double nNearInt; //nearest whole integer to n
        double jApprox; //approx value of j

        //n(*) = (J - J2000 - J0)
        n = (julianDay-J2K-JZERO) - longitude/360;
        nNearInt = Math.round(n);//nearest whole number
        jApprox = J2K + JZERO + longitude/360 + nNearInt;

        //re-calculate meanAnom for jApprox
        meanAnomDeg = meanAnomalyFromJulianDayTime(jApprox);

        double meanAnomRad = degreesToRadians(meanAnomDeg);
        double lSun = meanAnomDeg + ELOP + 180;
        //convert to range 0-360 if equal to or over 360. What about negative?
        if (lSun >= 360) {
            lSun = lSun - 360;
        }
        double lSunRad = degreesToRadians(lSun);

        //Jtransit J(*) + J1 sin M + J2 sin(2 Lsun)
        jTransit = jApprox + (JONE*(Math.sin(meanAnomRad))) + (JTWO*(Math.sin(2*lSunRad)));

        //repetition method to increase accuracy ~ jTransit = jTransit - hourAngleForJTransit/360
        do {
            System.out.println("                                             ");
            System.out.println("RECALCULATE HOUR ANGLE H (REPETITION METHOD)");
            meanAnomalyFromJulianDayTime(jTransit);
            equationCentre();
            eclipticalCoordinates();
            equatorialCoordinates();
            hourAngleFromJulianDayTime(jTransit);

            System.out.println("NEW HOUR ANGLE H -------------------------------------------------------------- " + hourAngle);// should reduce/tend to zero with iterations
            System.out.println("TIME OF TRANSIT (JULIAN DAY) -------------------------------------------------- " + jTransit);
            jTransit = jTransit - (hourAngle/360);
            System.out.println("TIME OF TRANSIT BY REPETITION METHOD ------------------------------------------ " + jTransit);
            //System.out.println("                                             ");
        } while (hourAngle > 0.0001);//close to zero
        return jTransit;
    }

    //calculate transit time (hours & minutes)
    public static void solarTransitTime() {
        System.out.println("");
        System.out.println("TIME OF TRANSIT (HOURS AND MINUTES) ------------------------------------------- " + julianDayToTime(jTransit, julianDayRaw));
        System.out.println("");
    }

    //calculate sunrise & sunset times
    public static void riseAndSet() {
        //constant for atmospheric refraction effect and apparent solar disc diameter (earth, degrees)
        final double H0 = -0.83;
        double correction = 0;

        //H = arccos((sin h??? ??? sin ?? sin ??)/(cos ?? cos ??))
        double H =  Math.acos((Math.sin(degreesToRadians(H0)) - (Math.sin(degreesToRadians(latitude)) *
                Math.sin(degreesToRadians(decSun))))/(Math.cos(degreesToRadians(latitude)) *
                Math.cos(degreesToRadians(decSun))));
        System.out.println("VALUE OF H (DEGREES) ---------------------------------------------------------- " + radiansToDegrees(H));

        //Jrise = jTransit - H/360 (*J3=1)
        jRise = jTransit - (radiansToDegrees(H)/360);//repetition method available     //first estimate
        System.out.println("TIME OF SUNRISE (JULIAN DAY) ---------FIRST ESTIMATE--------------------------- " + jRise);
        System.out.println("sunrise (HOURS AND MINUTES) first estimate ------------------------------------------- " + julianDayToTime(jRise, julianDayRaw));

        //start do while loop for repetition method here (sunrise)
        do {
            System.out.println("                   ");
            System.out.println("    rise             ");
            System.out.println("                   ");
            System.out.println("RECALCULATE H & DECSUN for JRISE (REPETITION METHOD)");
            meanAnomalyFromJulianDayTime(jRise);
            equationCentre();
            eclipticalCoordinates();
            equatorialCoordinates();
            //H = arccos((sin h??? ??? sin ?? sin ??)/(cos ?? cos ??))
            double hRise =  Math.acos((Math.sin(degreesToRadians(H0)) - (Math.sin(degreesToRadians(latitude)) *
                    Math.sin(degreesToRadians(decSun))))/(Math.cos(degreesToRadians(latitude)) *
                    Math.cos(degreesToRadians(decSun))));
            System.out.println("VALUE OF NEW HRISE (DEGREES) -------------------REPETITION METHOD-------------- " + radiansToDegrees(hRise));
            //Jrise = jTransit - H/360 (*J3=1)
            jRise = jTransit - (radiansToDegrees(hRise)/360);//repetition method     //repeat estimate
            System.out.println("TIME OF SUNRISE (JULIAN DAY) ---------REPEAT ESTIMATE-------------------------- " + jRise);
            correction = ((hRise - H)/360);
            System.out.println("CORRECTION = " + correction);
            System.out.println("                   ");
        } while (correction > 0.0001 || correction < -0.0001);

        //Jset = jTransit + H/360 (*J3=1)
        jSet = jTransit + (radiansToDegrees(H)/360);//repetition method available     //first estimate
        System.out.println("TIME OF SUNSET (JULIAN DAY) ---------------FIRST ESTIMATE---------------------- " + jSet);
        System.out.println("sunset (HOURS AND MINUTES) first estimate ------------------------------------------- " + julianDayToTime(jSet, julianDayRaw));

        //start do while loop for repetition method here (sunset)
        do {
            System.out.println("                   ");
            System.out.println("     set            ");
            System.out.println("                   ");
            System.out.println("RECALCULATE H & DECSUN for JSET (REPETITION METHOD)");
            meanAnomalyFromJulianDayTime(jSet);
            equationCentre();
            eclipticalCoordinates();
            equatorialCoordinates();
            //H = arccos((sin h??? ??? sin ?? sin ??)/(cos ?? cos ??))
            double hSet =  Math.acos((Math.sin(degreesToRadians(H0)) - (Math.sin(degreesToRadians(latitude)) *
                    Math.sin(degreesToRadians(decSun))))/(Math.cos(degreesToRadians(latitude)) *
                    Math.cos(degreesToRadians(decSun))));
            System.out.println("VALUE OF NEW HSET (DEGREES) -------------------REPETITION METHOD--------------- " + radiansToDegrees(hSet));
            //Jrise = jTransit - H/360 (*J3=1)
            jSet = jTransit + (radiansToDegrees(hSet)/360);//repetition method     //repeat estimate
            System.out.println("TIME OF SUNSET (JULIAN DAY) ---------REPEAT ESTIMATE--------------------------- " + jSet);
            correction = ((hSet - H)/360);
            System.out.println("CORRECTION = " + correction);
        } while (correction > 0.0001 || correction < -0.0001);
    }

    public static void riseAndSetTime() {
        System.out.println("  ");
        System.out.println("TIME OF SUNRISE (HOURS AND MINUTES) ------------------------------------------- " + julianDayToTime(jRise, julianDayRaw));
        System.out.println("TIME OF SUNSET (HOURS AND MINUTES) -------------------------------------------- " + julianDayToTime(jSet, julianDayRaw));
    }

    //convert Julian day number to time of day in hours and minutes
    public static String julianDayToTime (double jTime, double julianDayZero) {
        double timeHours = (jTime - julianDayZero)*24;
        double timeMinutes = (timeHours - (Math.floor(timeHours)));
        timeHours = timeHours - timeMinutes;
        timeMinutes = timeMinutes*60;
        if (timeMinutes >= 59.5) {
            timeMinutes = 0;
            timeHours = timeHours + 1;
        }
        if (timeHours < 0) {
            timeHours = timeHours + 24;
        }
        if (timeHours >= 24) {
            timeHours = timeHours - 24;
        }

        String outputString = dfw.format(timeHours) + "." + dfm.format(timeMinutes) + " (UTC)";
        return outputString;
    }

    //convert degrees to radians
    public static double degreesToRadians(double degrees)
    {
        double radians = degrees * 0.0174532925199;
        return radians;
    }

    //convert radians to degrees
    public static double radiansToDegrees(double radians)
    {
        double degrees = radians / 0.0174532925199;
        return degrees;
    }
}


