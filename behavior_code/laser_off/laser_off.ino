// Arduino Inhibition Pulse Code
// Written Lucius Wilmerding 5/7/26 - Grienberger Lab, Brandeis University

const int ledPin = 6;
const int aioPin = 5;

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

}