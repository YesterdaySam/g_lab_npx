// Arduino Inhibition Pulse Code
// Written Lucius Wilmerding 5/7/26 - Grienberger Lab, Brandeis University

const int ledPin = 6;
const int aioPin = 5;
const int rstPin = 2;
const int pulseWidth = 12;  //ms
const int pulseDelay = 25;  //ms total time between leading edge of each pulse (e.g. 30s up 30s down should be 60s pulseDelay)
const int nLapsPer = 5;
const int nCycles = 3;
const int rstinterval = 2000;  //ms
// int rstVal = 0;
int curLap = 0;
int curCycle = 1;
unsigned long currms = 0;
unsigned long prevms = 0;

void setup() {
  // put your setup code here, to run once:

  pinMode(ledPin, OUTPUT);
  pinMode(aioPin, OUTPUT);
  pinMode(LED_BUILTIN, OUTPUT);
  pinMode(rstPin, INPUT);
  digitalWrite(ledPin, LOW);
  digitalWrite(aioPin, LOW);
  digitalWrite(LED_BUILTIN, HIGH);
}

void loop() {
  // put your main code here, to run repeatedly:

  if (curCycle <= nCycles) {
    currms = millis();

    if (millis() - prevms >= rstinterval) {

      if (digitalRead(rstPin) == HIGH) {
        curLap++;
        prevms += rstinterval;
      }
    }

    if ((curLap > nLapsPer) && (curLap <= nLapsPer * 2)) {
      digitalWrite(LED_BUILTIN, LOW);
    }
    else if (curLap > nLapsPer * 2) {
      digitalWrite(LED_BUILTIN, HIGH);
      curLap = 0;
      curCycle++;
    }
  }
  else {
    digitalWrite(LED_BUILTIN, LOW);
  }
  //
  //while (currPulse < nPulses) {
  //
  //    delay(pulseDelay - pulseWidth);
  //
  //  digitalWrite(ledPin, HIGH);
  //  digitalWrite(aioPin, HIGH);
  //  digitalWrite(LED_BUILTIN, HIGH);
  //
  //   delay(pulseWidth);

  //    digitalWrite(ledPin, LOW);
  //    digitalWrite(aioPin, LOW);
  //    digitalWrite(LED_BUILTIN, LOW);

  //    currPulse++;
  //  }
}