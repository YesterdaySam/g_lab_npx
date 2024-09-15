// Arduino Phototag Code
// Written Lucius Wilmerding 9/15/2024 - Grienberger Lab, Brandeis University

const int ledPin = 6;
const int aioPin = 7;
const int pulseWidth = 15; //ms
const int pulseDelay = 1000; //ms every 5sec
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