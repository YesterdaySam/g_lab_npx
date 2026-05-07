// Arduino Inhibition Pulse Code
// Written Lucius Wilmerding 5/7/26 - Grienberger Lab, Brandeis University

const int ledPin = 6;
const int aioPin = 7;
const int pulseWidth = 30*1000; //ms
const int pulseDelay = 60*1000; //ms total time between leading edge of each pulse (e.g. 30s up 30s down should be 60s pulseDelay)
const int nPulses = 3; // 5 min
int currPulse = 0; 

void setup() {
  // put your setup code here, to run once:

  pinMode(ledPin, OUTPUT);
  pinMode(aioPin, OUTPUT);
  pinMode(LED_BUILTIN, OUTPUT);
  digitalWrite(ledPin, LOW);
  digitalWrite(aioPin, LOW);
  digitalWrite(LED_BUILTIN, LOW);
}

void loop() {
  // put your main code here, to run repeatedly:
  while (currPulse < nPulses) {

    delay(pulseDelay - pulseWidth);

    digitalWrite(ledPin, HIGH);
    digitalWrite(aioPin, HIGH);
    digitalWrite(LED_BUILTIN, HIGH);

    delay(pulseWidth);

    digitalWrite(ledPin, LOW);
    digitalWrite(aioPin, LOW);
    digitalWrite(LED_BUILTIN, LOW);

    currPulse++;
  }

}