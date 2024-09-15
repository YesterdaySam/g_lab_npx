// Arduino Flipper code to generate pseudo random timestamps
// From Tim Harris' Lab https://github.com/cortex-lab/neuropixels/wiki/Synchronization
// Modified by LKW; Grienberer lab, Brandeis University; 6/24/2024

const int poissonPin = 7;
const int minPoissonDur = 50; //ms
const int maxPoissonDur = 500; //ms

int currentPoissonState = 0;

void setup() {
  // put your setup code here, to run once:

  pinMode(poissonPin, OUTPUT);
  pinMode(LED_BUILTIN, OUTPUT);
  digitalWrite(poissonPin, LOW);
  currentPoissonState = LOW;
}

void loop() {
  // put your main code here, to run repeatedly:

  //code for flipper
  int poissonStateDur = random(minPoissonDur, maxPoissonDur);
  delay(poissonStateDur);
  if (currentPoissonState==LOW){
    currentPoissonState=HIGH;
  } 
  else {
    currentPoissonState=LOW;
  }
  digitalWrite(poissonPin, currentPoissonState);
  digitalWrite(LED_BUILTIN, currentPoissonState);
}