
const int ledPin = 5;
const int aioPin = 6;
const int poissonPin = 7;
const int minPoissonDur = 50; //ms
const int maxPoissonDur = 500; //ms
unsigned long tOpto;
unsigned long tFlip;
const int pulseWidth = 15; //ms
const int pulseDelay = 2000; //ms every 2sec
const int nPulses = 90; // 3 min
int currPulse = 0; 
int currentPoissonState = 0;

void setup() {
  // put your setup code here, to run once:
  pinMode(ledPin, OUTPUT);
  pinMode(aioPin, OUTPUT);
  pinMode(poissonPin, OUTPUT);
  pinMode(LED_BUILTIN, OUTPUT);
  
  digitalWrite(poissonPin, LOW);
  digitalWrite(ledPin, LOW);
  digitalWrite(aioPin, LOW);
  digitalWrite(LED_BUILTIN, LOW);
  
  currentPoissonState = LOW;
  tOpto = millis();
  tFlip = millis();
}

void loop() {
  // put your main code here, to run repeatedly
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
  
//  if (millis() - tFlip > poissonStateDur) {
//    tFlip = millis();
//    if (currentPoissonState == LOW) {
//      currentPoissonState = HIGH;
//    }
//    else {
//      currentPoissonState = LOW;
//    }
//    digitalWrite(poissonPin, currentPoissonState);
//  }
//
//  if (millis() - tFlip > poissonStateDur) {
//    if (currentPoissonState == LOW) {
//      currentPoissonState = HIGH;
//    }
//    else {
//      currentPoissonState = LOW;
//    }
//    digitalWrite(poissonPin, currentPoissonState);
//  }

  if (millis() - tOpto > pulseDelay & currPulse < nPulses) { 
    tOpto = millis();
    digitalWrite(ledPin, HIGH);
    digitalWrite(aioPin, HIGH);
    delay(pulseWidth);
    digitalWrite(ledPin, LOW);
    digitalWrite(aioPin, LOW);
    currPulse++;
  }
//
//  if (millis() - tOpto > pulseWidth) {
//    digitalWrite(ledPin, LOW);
//    //currPulse++;
//  }
}
