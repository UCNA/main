#include "ControlMenu.hh"
#include <cmath>

void calc_Add(std::deque<std::string>& args,std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(2,args,stack))
		return;
	float b = streamInteractor::popFloat(stack);
	float a = streamInteractor::popFloat(stack);
	printf("%g%+g = %g\n",a,b,a+b);
	stack.push(dtos(a+b));
}

void calc_Mul(std::deque<std::string>& args,std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(2,args,stack))
		return;
	float b = streamInteractor::popFloat(stack);
	float a = streamInteractor::popFloat(stack);
	printf("%g*%g = %g\n",a,b,a*b);
	stack.push(dtos(a*b));
}

void calc_Div(std::deque<std::string>& args,std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(2,args,stack))
		return;
	float b = streamInteractor::popFloat(stack);
	float a = streamInteractor::popFloat(stack);
	printf("%g/%g = %g\n",a,b,a/b);
	stack.push(dtos(a/b));
}

void calc_Pow(std::deque<std::string>& args,std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(2,args,stack))
		return;
	float b = streamInteractor::popFloat(stack);
	float a = streamInteractor::popFloat(stack);
	printf("%g^%g = %g\n",a,b,pow(a,b));
	stack.push(dtos(pow(a,b)));
}

void calc_Neg(std::deque<std::string>& args,std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(1,args,stack))
		return;
	float a = streamInteractor::popFloat(stack);
	printf("%g -> %g\n",a,-a);
	stack.push(dtos(-a));
}


void FPNCalc() {
	
	OptionsMenu CalcMain("FPN Calculator",false);

	inputRequester zAdd("Add x+y",&calc_Add);
	zAdd.addArg("x","",&CalcMain);
	zAdd.addArg("y","",&CalcMain);
	inputRequester zMul("Multiply x*y",&calc_Mul);
	zMul.addArg("x","",&CalcMain);
	zMul.addArg("y","",&CalcMain);
	inputRequester zDiv("Divide x/y",&calc_Div);
	zDiv.addArg("x","",&CalcMain);
	zDiv.addArg("y","",&CalcMain);
	inputRequester zPow("Exponent x^y",&calc_Pow);
	zPow.addArg("x","",&CalcMain);
	zPow.addArg("y","",&CalcMain);
	inputRequester zNeg("Negate -x",&calc_Neg);
	zNeg.addArg("x","",&CalcMain);
	
	streamInteractor identity;
	
	CalcMain.addChoice(&zAdd,"+");
	CalcMain.addChoice(&zMul,"*");
	CalcMain.addChoice(&zDiv,"/");
	CalcMain.addChoice(&zPow,"^");
	CalcMain.addChoice(&zNeg,"-");
	
	CalcMain.setCatchall(&identity);
	
	std::deque<std::string> d;
	std::stack<std::string> s;
	while(!d.size() || !startsWith(d.front(),NameSelector::exit_control)) {
		CalcMain.doIt(d,s);
		menutils_PrintStack(d,s);
	}
}

void RPNCalc() {
	
	OptionsMenu CalcMain("RPN Calculator",false);
	
	inputRequester zAdd("Add x+y",&calc_Add);
	inputRequester zMul("Multiply x*y",&calc_Mul);
	inputRequester zDiv("Divide x/y",&calc_Div);
	inputRequester zPow("Exponent x^y",&calc_Pow);
	inputRequester zNeg("Negate -x",&calc_Neg);
	
	inputRequester zDrop("Drop top item",&menutils_Drop);
	inputRequester zDup("Duplicate top item",&menutils_Dup);
	inputRequester zSwap("Swap top 2 items",&menutils_Swap);
	inputRequester zRot("Rotate top n items",&menutils_Rot);
	inputRequester zClear("Clear stack",&menutils_ClearStack);
	inputRequester zExec("Execute commands",&menutils_Exec);
	inputRequester zSelect("Select c?a:b",&menutils_Select);
	
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
	
	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	CalcMain.addChoice(&exitMenu,"x");
	CalcMain.addSynonym("x","quit");
	CalcMain.addSynonym("x","bye");
	
	streamInteractor identity;
	CalcMain.setCatchall(&identity);
	
	std::deque<std::string> d;
	std::stack<std::string> s;
	while(!d.size() || !startsWith(d.front(),NameSelector::exit_control) && !startsWith(d.front(),NameSelector::barf_control)) {
		CalcMain.doIt(d,s);
		menutils_PrintStack(d,s);
		printf("\n");
	}
	menutils_PrintQue(d,s);
}

int main(int argc, char *argv[]) {
	printf("Starting calculator...\n");
	RPNCalc();
	printf("Goodbye.\n");
	return 0;
}
