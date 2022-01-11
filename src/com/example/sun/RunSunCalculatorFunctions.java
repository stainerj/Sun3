package com.example.sun;

public class RunSunCalculatorFunctions {

    public static void main(String[] args) {

        SunCalculatorFunctions suncalc = new SunCalculatorFunctions();
        suncalc.setTimeZone(0);
        suncalc.dateTimeToJulianDay(2021, 12, 12, 8.00);
        suncalc.setLatLong(54.6, -5.9);
        suncalc.meanAnomaly();
        suncalc.equationCentre();
        suncalc.eclipticalCoordinates();
        suncalc.equatorialCoordinates();
        suncalc.observer();
        suncalc.solarTransit();
        suncalc.solarTransitTime();
        suncalc.riseAndSet();
        suncalc.riseAndSetTime();

        //getters?
        //return object of all output data?
    }
}
