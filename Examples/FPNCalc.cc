#include "ControlMenu.hh"
#include <cmath>

void calc_Add(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(2))
		return;
	float b = S->popFloat();
	float a = S->popFloat();
	printf("%g%+g = %g\n",a,b,a+b);
	S->mystack->push(dtos(a+b));
}

void calc_Mul(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(2))
		return;
	float b = S->popFloat();
	float a = S->popFloat();
	printf("%g*%g = %g\n",a,b,a*b);
	S->mystack->push(dtos(a*b));
}

void calc_Div(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(2))
		return;
	float b = S->popFloat();
	float a = S->popFloat();
	printf("%g/%g = %g\n",a,b,a/b);
	S->mystack->push(dtos(a/b));
}

void calc_Pow(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(2))
		return;
	float b = S->popFloat();
	float a = S->popFloat();
	printf("%g^%g = %g\n",a,b,pow(a,b));
	S->mystack->push(dtos(pow(a,b)));
}

void calc_Neg(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(1))
		return;
	float a = S->popFloat();
	printf("%g -> %g\n",a,-a);
	S->mystack->push(dtos(-a));
}


void FPNCalc() {
	
	OptionsMenu CalcMain("FPN Calculator",false);

	InputRequester zAdd("Add x+y",&calc_Add);
	zAdd.addArg("x","","",&CalcMain);
	zAdd.addArg("y","","",&CalcMain);
	InputRequester zMul("Multiply x*y",&calc_Mul);
	zMul.addArg("x","","",&CalcMain);
	zMul.addArg("y","","",&CalcMain);
	InputRequester zDiv("Divide x/y",&calc_Div);
	zDiv.addArg("x","","",&CalcMain);
	zDiv.addArg("y","","",&CalcMain);
	InputRequester zPow("Exponent x^y",&calc_Pow);
	zPow.addArg("x","","",&CalcMain);
	zPow.addArg("y","","",&CalcMain);
	InputRequester zNeg("Negate -x",&calc_Neg);
	zNeg.addArg("x","","",&CalcMain);
	
	StreamInteractor identity;
	
	CalcMain.addChoice(&zAdd,"+");
	CalcMain.addChoice(&zMul,"*");
	CalcMain.addChoice(&zDiv,"/");
	CalcMain.addChoice(&zPow,"^");
	CalcMain.addChoice(&zNeg,"-");
	
	CalcMain.setCatchall(&identity);
	
	std::deque<std::string> d;
	std::stack<std::string> s;
	CalcMain.mydeque = &d;
	CalcMain.mystack = &s;
	while(!d.size() || !startsWith(d.front(),NameSelector::exit_control)) {
		CalcMain.doIt();
		menutils_PrintStack(&CalcMain);
	}
}

void RPNCalc() {
	
	OptionsMenu CalcMain("RPN Calculator",false);
	
	InputRequester zAdd("Add x+y",&calc_Add);
	InputRequester zMul("Multiply x*y",&calc_Mul);
	InputRequester zDiv("Divide x/y",&calc_Div);
	InputRequester zPow("Exponent x^y",&calc_Pow);
	InputRequester zNeg("Negate -x",&calc_Neg);
	
	InputRequester zDrop("Drop top item",&menutils_Drop);
	InputRequester zDup("Duplicate top item",&menutils_Dup);
	InputRequester zSwap("Swap top 2 items",&menutils_Swap);
	InputRequester zRot("Rotate top n items",&menutils_Rot);
	InputRequester zClear("Clear stack",&menutils_ClearStack);
	InputRequester zExec("Execute commands",&menutils_Exec);
	InputRequester zSelect("Select c?a:b",&menutils_Select);
	
	CalcMain.addChoice(&zAdd,"+");
	CalcMain.addChoice(&zMul,"*");
	CalcMain.addChoice(&zDiv,"/");
	CalcMain.addChoice(&zPow,"^");
	CalcMain.addChoice(&zNeg,"-");
	
	CalcMain.addChoice(&zDrop,"drop");
	CalcMain.addChoice(&zDup,"dup");
	CalcMain.addChoice(&zSwap,"swap");
	CalcMain.addChoice(&zRot,"rot");
	CalcMain.addChoice(&zClear,"clear");
	CalcMain.addChoice(&zExec,"exec");
	CalcMain.addChoice(&zSelect,"sel");
	
	InputRequester exitMenu("Exit Menu",&menutils_Exit);
	CalcMain.addChoice(&exitMenu,"x");
	CalcMain.addSynonym("x","quit");
	CalcMain.addSynonym("x","bye");
	
	StreamInteractor identity;
	CalcMain.setCatchall(&identity);
	
	std::deque<std::string> d;
	std::stack<std::string> s;
	CalcMain.mydeque = &d;
	CalcMain.mystack = &s;
	while(!d.size() || (!startsWith(d.front(),NameSelector::exit_control) && !startsWith(d.front(),NameSelector::barf_control))) {
		CalcMain.doIt();
		menutils_PrintStack(&CalcMain);
		printf("\n");
	}
	menutils_PrintQue(&CalcMain);
}

int main(int, char **) {
	printf("Starting calculator...\n");
	RPNCalc();
	printf("Goodbye.\n");
	return 0;
}
