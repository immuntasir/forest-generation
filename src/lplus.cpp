#include "lplus.h"
#include <fstream>
#include "R3Mesh.h"
using namespace std;
LPlusSystem::LPlusSystem(R3Mesh * m)
:LSystem(m),isPlus(false)
{
}
double leaf_size;
string LPlusSystem::generateFromFile(const char * filename,const int iterationsOverride, double leafSize)
{
	int l=strlen(filename);
	leaf_size = leafSize;
	if (strcmp(filename+l-3,"l++")==0) //this is an l++
	{
		isPlus=true;
	}
	return LSystem::generateFromFile(filename,iterationsOverride);

}
void LPlusSystem::setTurtlePosition(float x, float y, float z){
	turtle.position = R3Vector(x, 0, z);
	turtle.direction = R3Vector(0, 1, 0);
}

float rand_gen() {
   // return a uniformly distributed random value
   return ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
}
float normalRandom(float mu, float sigma) {
   // return a normally distributed random value
   float v1=rand_gen();
   float v2=rand_gen();
   float seed = cos(2*3.14*v2)*sqrt(-2.*log(v1));
   if (seed>2) seed = 2;
   if (seed<-2) seed = -2;
   return sigma * seed + mu;
}

void LPlusSystem::run(const char command,const float param)
{
	if (!isPlus)
		return LSystem::run(command,param);

	float co=defaultCoefficient;
	float num=param;
	if (num==1)
		num*=co;
	float randomNumber;
	switch (command)
	{
		case '+':
		turtle.turnLeft(num);
		break;
		case '-':
		turtle.turnRight(num);
		break;
		case '&':
		randomNumber = normalRandom(num, num*.1);
		turtle.pitchDown(randomNumber);
		break;
		case '^':
		randomNumber = normalRandom(num, num*.1);
		turtle.pitchUp(randomNumber);
		break;
		case '<':
		turtle.thicken(num);
		break;
		case '\\':
		randomNumber = normalRandom(num, num*.1);
		turtle.rollLeft(randomNumber);
		break;
		case '/':
		////////////////////////////////////////  L++ ////////////////////////////////////////////////
		randomNumber = normalRandom(num, num*.1);
		turtle.rollRight(randomNumber);
		break;
		case '>':
		randomNumber = normalRandom(num, num*.1);
		turtle.narrow(randomNumber);
		break;
		case '%':
		turtle.setReduction(param);
		break;
		case '=':
		turtle.setThickness(param);
		break;
		case '|':
		turtle.turn180(param);
		break;
		case '*':
		randomNumber = normalRandom(leaf_size, leaf_size*.15);
		turtle.drawLeaf(randomNumber);
		break;
		case 'F':
		case 'f':
		randomNumber = normalRandom(param, param*.2);
		turtle.draw(param);
		case 'G':
		case 'g':
		turtle.move(param);
		break;
		case '[':
			turtle.save();
			break;
		case ']':
			turtle.restore();
			break;
		default:
		;
	}

}
