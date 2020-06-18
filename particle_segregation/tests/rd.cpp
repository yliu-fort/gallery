#include <iostream>
#include <random>

using namespace std;

int main()
{
    // Matlab seed=5418, 'twister'
    //0.828444499439415
    //0.104379647114743
    //0.342207534058427

    float val;
    std::mt19937 rng(77777);
    std::uniform_real_distribution<float> urd(0, 1);
    val = urd(rng);
    cout << val << endl;
    /* skip alternate values */ rng();
    val = urd(rng);
    cout << val << endl;
    /* skip alternate values */ rng();
    val = urd(rng);
    cout << val << endl;
}
