#ifndef COLORMAP_H
#define COLORMAP_H

#include <iostream>
#include <vector>

//GLEW
#define GLEW_STATIC
#include <GL/glew.h>

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

// initialize color vector and bind to gltexture1d
namespace Colormap {

// handle in cpu
static std::vector<glm::vec3> _colorMap;
// handle in gpu
static unsigned int *_colorMapHandle = nullptr;
// Texture binding point
static int _id = 10;

// _func() is undesired to be called from out of the scope
// destroy method
void _release()
{
    glDeleteTextures(1, _colorMapHandle);
    delete _colorMapHandle;
    _colorMapHandle = nullptr;
}

// APIs: expected to be called from outerscope
// initialization method
void GLinit()
{
    if(_colorMapHandle) _release();

    // Allocate and bind texture buffer
    _colorMapHandle = new unsigned int;
    glGenTextures(1, _colorMapHandle);
    glBindTexture(GL_TEXTURE_1D, *_colorMapHandle);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB16, int(_colorMap.size()), 0, GL_RGB, GL_FLOAT, &_colorMap[0]);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
}

void Bind()
{
    // Check if colormap handle is configured
    if(!_colorMapHandle) return;

    // Bind to texture
    glActiveTexture(GL_TEXTURE0+_id);
    glBindTexture(GL_TEXTURE_1D, *_colorMapHandle);
}

void Bind(int id)
{
    // Check if colormap handle is configured
    if(!_colorMapHandle) return;

    // Bind to texture
    glActiveTexture(GL_TEXTURE0+id);
    glBindTexture(GL_TEXTURE_1D, *_colorMapHandle);
}

int ID()
{
    return _id;
}

// switch methods
void Rainbow()
{
    _colorMap.clear();

    // assign colormap
    unsigned int mapResolution = 64;
    for(unsigned int i = 0; i < mapResolution; i++)
    {
        /*plot short rainbow RGB*/
        double f = i/(mapResolution - 1.0);
        double a=(1-f)/0.25;    //invert and group
        int X=int(floor(a));    //this is the integer part
        int Y=int(floor(255*(a-X))); //fractional part from 0 to 255
        float r=0.0f,g=0.0f,b=0.0f;
        switch(X)
        {
        case 0: r=255;g=Y;b=0;break; // RED f = 1
        case 1: r=255-Y;g=255;b=0;break;
        case 2: r=0;g=255;b=Y;break;
        case 3: r=0;g=255-Y;b=255;break;
        case 4: r=0;g=0;b=255;break; // BLUE f = 0
        }

        glm::vec3 color(r/255.0f, g/255.0f, b/255.0f);
        _colorMap.push_back(color);
    }

    GLinit();
}

void Cool_warm()
{
    _colorMap.clear();

    /*plot short rainbow RGB*/
    _colorMap.push_back(glm::vec3(0.3347908501	,0.2830843703	,0.7564952199));
    _colorMap.push_back(glm::vec3(0.3518280609	,0.3116758913	,0.7812772411));
    _colorMap.push_back(glm::vec3(0.3688881862	,0.3399620141	,0.8050035716));
    _colorMap.push_back(glm::vec3(0.3859900446	,0.3679393904	,0.8276127105));
    _colorMap.push_back(glm::vec3(0.4031497486	,0.3955898511	,0.8490463938));
    _colorMap.push_back(glm::vec3(0.4203805022	,0.422885173     ,0.8692497346));
    _colorMap.push_back(glm::vec3(0.437692408   ,0.4497901767	,0.888171355 ));
    _colorMap.push_back(glm::vec3(0.4550922946	,0.4762648209	,0.9057635071));
    _colorMap.push_back(glm::vec3(0.472583575   ,0.5022656703	,0.9219821839));
    _colorMap.push_back(glm::vec3(0.4901661419	,0.5277469573	,0.9367872181));
    _colorMap.push_back(glm::vec3(0.5078363019	,0.552661374     ,0.95014237  ));
    _colorMap.push_back(glm::vec3(0.5255867527	,0.5769606806	,0.9620154023));
    _colorMap.push_back(glm::vec3(0.5434066     ,0.6005961848	,0.972378143 ));
    _colorMap.push_back(glm::vec3(0.5612814146	,0.6235191302	,0.9812065351));
    _colorMap.push_back(glm::vec3(0.5791933257	,0.6456810176	,0.9884806737));
    _colorMap.push_back(glm::vec3(0.597121147   ,0.6670338772	,0.9941848299));
    _colorMap.push_back(glm::vec3(0.6150405306	,0.6875305036	,0.9983074607));
    _colorMap.push_back(glm::vec3(0.6329241457	,0.707124662     ,1.0008412073));
    _colorMap.push_back(glm::vec3(0.6507418762	,0.725771272     ,1.0017828782));
    _colorMap.push_back(glm::vec3(0.6684610336	,0.7434265724	,1.0011334207));
    _colorMap.push_back(glm::vec3(0.6860465797	,0.760048273     ,0.9988978787));
    _colorMap.push_back(glm::vec3(0.7034613582	,0.7755956914	,0.9950853369));
    _colorMap.push_back(glm::vec3(0.7206663277	,0.7900298814	,0.9897088538));
    _colorMap.push_back(glm::vec3(0.7376207974	,0.80331375      ,0.9827853805));
    _colorMap.push_back(glm::vec3(0.7542826598	,0.8154121667	,0.9743356681));
    _colorMap.push_back(glm::vec3(0.7706086203	,0.8262920645	,0.9643841623));
    _colorMap.push_back(glm::vec3(0.7865544213	,0.8359225335	,0.9529588871));
    _colorMap.push_back(glm::vec3(0.8020750595	,0.8442749067	,0.9400913151));
    _colorMap.push_back(glm::vec3(0.8171249966	,0.8513228389	,0.9258162279));
    _colorMap.push_back(glm::vec3(0.8316583613	,0.8570423782	,0.9101715642));
    _colorMap.push_back(glm::vec3(0.8456291429	,0.8614120303	,0.8931982563));
    _colorMap.push_back(glm::vec3(0.8589913773	,0.8644128151	,0.8749400557));
    _colorMap.push_back(glm::vec3(0.8729974744	,0.8618078828	,0.8542093518));
    _colorMap.push_back(glm::vec3(0.8871428255	,0.8536182997	,0.8315274631));
    _colorMap.push_back(glm::vec3(0.8999602396	,0.8440896036	,0.8083454206));
    _colorMap.push_back(glm::vec3(0.9114398012	,0.8332425531	,0.7847243245));
    _colorMap.push_back(glm::vec3(0.9215721144	,0.821100097     ,0.7607246682));
    _colorMap.push_back(glm::vec3(0.9303485868	,0.8076872307	,0.7364061852));
    _colorMap.push_back(glm::vec3(0.9377616713	,0.7930308341	,0.7118277037));
    _colorMap.push_back(glm::vec3(0.9438050729	,0.7771594886	,0.6870470086));
    _colorMap.push_back(glm::vec3(0.9484739274	,0.7601032682	,0.6621207123));
    _colorMap.push_back(glm::vec3(0.9517649556	,0.7418935004	,0.637104134 ));
    _colorMap.push_back(glm::vec3(0.9536765968	,0.7225624877	,0.6120511882));
    _colorMap.push_back(glm::vec3(0.954209125   ,0.7021431827	,0.5870142815));
    _colorMap.push_back(glm::vec3(0.953364749   ,0.6806688005	,0.56204422  ));
    _colorMap.push_back(glm::vec3(0.9511476991	,0.6581723531	,0.5371901261));
    _colorMap.push_back(glm::vec3(0.9475643007	,0.6346860784	,0.5124993657));
    _colorMap.push_back(glm::vec3(0.9426230369	,0.6102407295	,0.4880174858));
    _colorMap.push_back(glm::vec3(0.9363346003	,0.5848646713	,0.4637881639));
    _colorMap.push_back(glm::vec3(0.9287119347	,0.5585827109	,0.4398531685));
    _colorMap.push_back(glm::vec3(0.9197702679	,0.5314145461	,0.4162523326));
    _colorMap.push_back(glm::vec3(0.9095271348	,0.5033726591	,0.3930235404));
    _colorMap.push_back(glm::vec3(0.8980023922	,0.4744593786	,0.3702027279));
    _colorMap.push_back(glm::vec3(0.8852182252	,0.4446626581	,0.3478239008));
    _colorMap.push_back(glm::vec3(0.8711991452	,0.4139498021	,0.3259191675));
    _colorMap.push_back(glm::vec3(0.8559719797	,0.3822577695	,0.3045187937));
    _colorMap.push_back(glm::vec3(0.839565854   ,0.3494774788	,0.2836512771));
    _colorMap.push_back(glm::vec3(0.8220121659	,0.3154269254	,0.2633434476));
    _colorMap.push_back(glm::vec3(0.8033445511	,0.2798017272	,0.2436205956));
    _colorMap.push_back(glm::vec3(0.7835988421	,0.2420751664	,0.2245066322));
    _colorMap.push_back(glm::vec3(0.7628130186	,0.2012677301	,0.2060242863));
    _colorMap.push_back(glm::vec3(0.7410271511	,0.1552973138	,0.1881953425));
    _colorMap.push_back(glm::vec3(0.7182833362	,0.0983399315	,0.1710409249));
    _colorMap.push_back(glm::vec3(0.6946256248	,0.0029646105	,0.1545818283));

    GLinit();
}

void Cool_warm_bent()
{
    _colorMap.clear();

    /*plot short rainbow RGB*/
    _colorMap.push_back(glm::vec3(0.334790850135,	0.283084370265	,0.756495219864));
    _colorMap.push_back(glm::vec3(0.346186237067,	0.306343359743	,0.766844849485));
    _colorMap.push_back(glm::vec3(0.358070886588,	0.329227925965	,0.77681840476 ));
    _colorMap.push_back(glm::vec3(0.370455145382,	0.35181470768	,0.786423420089));
    _colorMap.push_back(glm::vec3(0.383348014038,	0.37416061155	,0.795667861721));
    _colorMap.push_back(glm::vec3(0.396757258725,	0.396308591789	,0.80456011256 ));
    _colorMap.push_back(glm::vec3(0.410689526369,	0.41829146389	,0.81310895801 ));
    _colorMap.push_back(glm::vec3(0.425150459774,	0.440134504348	,0.821323572852));
    _colorMap.push_back(glm::vec3(0.44014480975	,   0.461857272238	,0.829213509087));
    _colorMap.push_back(glm::vec3(0.455676541901,	0.483474916026	,0.83678868475 ));
    _colorMap.push_back(glm::vec3(0.471748936333,	0.504999130445	,0.84405937364 ));
    _colorMap.push_back(glm::vec3(0.488364679077,	0.526438869821	,0.85103619598 ));
    _colorMap.push_back(glm::vec3(0.505525944528,	0.547800888339	,0.857730109955));
    _colorMap.push_back(glm::vec3(0.523234468569,	0.569090155062	,0.864152404129));
    _colorMap.push_back(glm::vec3(0.541491612419,	0.590310176852	,0.870314690716));
    _colorMap.push_back(glm::vec3(0.560298417474,	0.611463252526	,0.876228899672));
    _colorMap.push_back(glm::vec3(0.579655651623,	0.632550675052	,0.881907273604));
    _colorMap.push_back(glm::vec3(0.599563847625,	0.653572893969	,0.887362363441));
    _colorMap.push_back(glm::vec3(0.620023334244,	0.674529647069	,0.892607024857));
    _colorMap.push_back(glm::vec3(0.641034260857,	0.695420068071	,0.897654415393));
    _colorMap.push_back(glm::vec3(0.662596616263,	0.716242775394	,0.902517992236));
    _colorMap.push_back(glm::vec3(0.68471024242	,   0.736995945924	,0.907211510607));
    _colorMap.push_back(glm::vec3(0.707374843784,	0.757677376779	,0.911749022691));
    _colorMap.push_back(glm::vec3(0.730589992886,	0.77828453743	,0.916144877045));
    _colorMap.push_back(glm::vec3(0.754355132745,	0.798814613992	,0.920413718408));
    _colorMap.push_back(glm::vec3(0.778669576636,	0.819264547168	,0.924570487816));
    _colorMap.push_back(glm::vec3(0.803532505694,	0.839631064994	,0.928630422934));
    _colorMap.push_back(glm::vec3(0.828942964771,	0.859910711315	,0.932609058496));
    _colorMap.push_back(glm::vec3(0.85489985692	,   0.880099870754	,0.936522226723));
    _colorMap.push_back(glm::vec3(0.881401936818,	0.900194790782	,0.940386057604));
    _colorMap.push_back(glm::vec3(0.908447803397,	0.920191601388	,0.944216978887));
    _colorMap.push_back(glm::vec3(0.93603589193	,   0.940086332774	,0.948031715644));
    _colorMap.push_back(glm::vec3(0.946615483498,	0.938119667501	,0.932368364096));
    _colorMap.push_back(glm::vec3(0.939831686063,	0.914317074266	,0.897622658829));
    _colorMap.push_back(glm::vec3(0.933096113635,	0.890443008017	,0.863413341256));
    _colorMap.push_back(glm::vec3(0.926386773449,	0.866497010372	,0.82975302552 ));
    _colorMap.push_back(glm::vec3(0.919683081312,	0.8424780094	,0.796653973917));
    _colorMap.push_back(glm::vec3(0.912965817319,	0.81838423549	,0.764128104993));
    _colorMap.push_back(glm::vec3(0.906217075288,	0.794213120897	,0.732187004021));
    _colorMap.push_back(glm::vec3(0.899420207204,	0.769961179056	,0.700841936132));
    _colorMap.push_back(glm::vec3(0.892559763886,	0.745623858639	,0.67010386245 ));
    _colorMap.push_back(glm::vec3(0.885621432938,	0.721195365814	,0.639983459634));
    _colorMap.push_back(glm::vec3(0.878591974951,	0.696668446156	,0.610491143264));
    _colorMap.push_back(glm::vec3(0.871459158782,	0.67203411482	,0.5816370956  ));
    _colorMap.push_back(glm::vec3(0.864211696663,	0.647281319717	,0.553431298295));
    _colorMap.push_back(glm::vec3(0.856839179744,	0.6223965169	,0.525883570739));
    _colorMap.push_back(glm::vec3(0.849332014629,	0.597363129454	,0.499003614757));
    _colorMap.push_back(glm::vec3(0.841681361344,	0.57216084958	,0.472801066498));
    _colorMap.push_back(glm::vec3(0.833879073111,	0.5467647263	,0.447285556386));
    _colorMap.push_back(glm::vec3(0.825917638246,	0.521143954835	,0.42246677804 ));
    _colorMap.push_back(glm::vec3(0.817790124418,	0.495260242538	,0.398354567065));
    _colorMap.push_back(glm::vec3(0.809490125473,	0.469065560061	,0.3749589905  ));
    _colorMap.push_back(glm::vec3(0.801011710969,	0.442498976667	,0.352290447448));
    _colorMap.push_back(glm::vec3(0.79234937855	,   0.415482089745	,0.330359780961));
    _colorMap.push_back(glm::vec3(0.783498009224,	0.387912219692	,0.309178400404));
    _colorMap.push_back(glm::vec3(0.77445282562	,   0.359651901974	,0.288758412186));
    _colorMap.push_back(glm::vec3(0.765209353253,	0.330511927783	,0.269112754642));
    _colorMap.push_back(glm::vec3(0.755763384824,	0.30022242597	,0.250255329661));
    _colorMap.push_back(glm::vec3(0.746110947561,	0.268379963633	,0.23220111898 ));
    _colorMap.push_back(glm::vec3(0.736248273613,	0.234341304956	,0.214966266402));
    _colorMap.push_back(glm::vec3(0.726171773492,	0.196980137192	,0.198568098098));
    _colorMap.push_back(glm::vec3(0.715878012582,	0.154005981698	,0.183025041292));
    _colorMap.push_back(glm::vec3(0.705363690711,	0.0992176945762	,0.168356387181));
    _colorMap.push_back(glm::vec3(0.694625624821,0.00296461045768	,0.154581828278));

    GLinit();
}

void Viridis()
{
    _colorMap.clear();

    /*plot short rainbow RGB*/
    _colorMap.push_back(glm::vec3(0.282364552929017	,0                   ,0.331010194011805));
    _colorMap.push_back(glm::vec3(0.290264610268739	,0.00290166842697	,0.355012546339847));
    _colorMap.push_back(glm::vec3(0.29702948940838	,0.026329692233541	,0.377924541755913));
    _colorMap.push_back(glm::vec3(0.302631925424469	,0.053862907648956	,0.399626263232077));
    _colorMap.push_back(glm::vec3(0.307039918387932	,0.078713297374715	,0.420005358571776));
    _colorMap.push_back(glm::vec3(0.310233463566113	,0.101727325918509	,0.438954588437629));
    _colorMap.push_back(glm::vec3(0.312201373981498	,0.123749572612164	,0.456382208465964));
    _colorMap.push_back(glm::vec3(0.312951437700377	,0.145157350614433	,0.472222572547925));
    _colorMap.push_back(glm::vec3(0.312512607615915	,0.16612316070645	,0.486434918529507));
    _colorMap.push_back(glm::vec3(0.31093764455559	,0.186719832950713	,0.499013125962462));
    _colorMap.push_back(glm::vec3(0.308303989359403	,0.206966530927814	,0.509988416696046));
    _colorMap.push_back(glm::vec3(0.304710416967093	,0.226855800013103	,0.519427465222488));
    _colorMap.push_back(glm::vec3(0.300275109653277	,0.246369197021837	,0.527427552426432));
    _colorMap.push_back(glm::vec3(0.295131439828935	,0.265485180439366	,0.534112823026301));
    _colorMap.push_back(glm::vec3(0.289414399393372	,0.284190100450387	,0.539622289939683));
    _colorMap.push_back(glm::vec3(0.283259897633911	,0.302476118715382	,0.544103261023281));
    _colorMap.push_back(glm::vec3(0.276799621009854	,0.320345940507622	,0.5477066116626  ));
    _colorMap.push_back(glm::vec3(0.270145565408322	,0.337814625265202	,0.550571633159153));
    _colorMap.push_back(glm::vec3(0.263395095437236	,0.354902174182077	,0.552826728927636));
    _colorMap.push_back(glm::vec3(0.256627021197388	,0.371637422456423	,0.554587351226669));
    _colorMap.push_back(glm::vec3(0.249889684649734	,0.38805519013154	,0.555944447227502));
    _colorMap.push_back(glm::vec3(0.243211275057896	,0.404191590516798	,0.556972229077804));
    _colorMap.push_back(glm::vec3(0.236596479528174	,0.420087300026636	,0.557723111447191));
    _colorMap.push_back(glm::vec3(0.230018460633008	,0.435782119037193	,0.558229794787176));
    _colorMap.push_back(glm::vec3(0.223435922319333	,0.451315765364287	,0.558505056412238));
    _colorMap.push_back(glm::vec3(0.216792930821148	,0.466723486929379	,0.558544840752747));
    _colorMap.push_back(glm::vec3(0.210017477313124	,0.482038854278792	,0.558326659896291));
    _colorMap.push_back(glm::vec3(0.203033654502342	,0.497293987417666	,0.557811096326849));
    _colorMap.push_back(glm::vec3(0.195775301016666	,0.512512430604381	,0.556950100370821));
    _colorMap.push_back(glm::vec3(0.188187099057271	,0.527712648226893	,0.555681981825937));
    _colorMap.push_back(glm::vec3(0.180236959840526	,0.542910294909534	,0.553936330846346));
    _colorMap.push_back(glm::vec3(0.171944648231779	,0.558111888569955	,0.551642155173208));
    _colorMap.push_back(glm::vec3(0.163379666628041	,0.573320182095253	,0.548716977566472));
    _colorMap.push_back(glm::vec3(0.154712976838942	,0.58852981438377	,0.545084527878414));
    _colorMap.push_back(glm::vec3(0.146228489262813	,0.603733399040755	,0.540664156524644));
    _colorMap.push_back(glm::vec3(0.138382992372823	,0.618912709364973	,0.535382365747223));
    _colorMap.push_back(glm::vec3(0.1318268679528     ,0.634050951144566	,0.529166506376451));
    _colorMap.push_back(glm::vec3(0.127421964028518	,0.649122125180978	,0.521950134762619));
    _colorMap.push_back(glm::vec3(0.126162847462407	,0.664098432326387	,0.513672984321242));
    _colorMap.push_back(glm::vec3(0.129002817951699	,0.678948620438736	,0.504279884061317));
    _colorMap.push_back(glm::vec3(0.136601245268402	,0.693639426460058	,0.493722788967818));
    _colorMap.push_back(glm::vec3(0.149154621822552	,0.708133352836865	,0.481958321528951));
    _colorMap.push_back(glm::vec3(0.166406002460385	,0.722391786823486	,0.468951724515542));
    _colorMap.push_back(glm::vec3(0.187854550610124	,0.736373835143009	,0.454666444833591));
    _colorMap.push_back(glm::vec3(0.212871474712786	,0.750039786447896	,0.439079529327591));
    _colorMap.push_back(glm::vec3(0.240901722041613	,0.76334487740616	,0.422163555778038));
    _colorMap.push_back(glm::vec3(0.271480477780217	,0.776246820025822	,0.403905416079382));
    _colorMap.push_back(glm::vec3(0.304247328808078	,0.788702440746788	,0.384278183697825));
    _colorMap.push_back(glm::vec3(0.338911907701283	,0.80066809511891	,0.36327062820589 ));
    _colorMap.push_back(glm::vec3(0.375246625025813	,0.812101199712745	,0.340864183867427));
    _colorMap.push_back(glm::vec3(0.413064476351507	,0.822962348627194	,0.317041408353394));
    _colorMap.push_back(glm::vec3(0.452190206057646	,0.833214002282859	,0.291792440200484));
    _colorMap.push_back(glm::vec3(0.492471275324835	,0.842826009221976	,0.265094198612632));
    _colorMap.push_back(glm::vec3(0.533740231244458	,0.851774167176046	,0.236941433193716));
    _colorMap.push_back(glm::vec3(0.575818711351207	,0.860049030662412	,0.207338402725059));
    _colorMap.push_back(glm::vec3(0.618512934051516	,0.867656904847372	,0.176318808263156));
    _colorMap.push_back(glm::vec3(0.661598108005148	,0.87462285455815	,0.144012124901709));
    _colorMap.push_back(glm::vec3(0.704830690454684	,0.880997166055619	,0.110791626635145));
    _colorMap.push_back(glm::vec3(0.747939656746893	,0.886858280605159	,0.077798082377023));
    _colorMap.push_back(glm::vec3(0.790657955615173	,0.892306987639634	,0.048645611984045));
    _colorMap.push_back(glm::vec3(0.832724255401618	,0.89746720147354	,0.033814584065587));
    _colorMap.push_back(glm::vec3(0.873918125624655	,0.902469498021461	,0.04447574610246 ));
    _colorMap.push_back(glm::vec3(0.914039449780679	,0.907456749130336	,0.074152638998232));
    _colorMap.push_back(glm::vec3(0.952999453291615	,0.91254523282901	,0.110858769093613));

    GLinit();
}

void Plasma()
{
    _colorMap.clear();

    /*plot short rainbow RGB*/
    _colorMap.push_back(glm::vec3(0.185001262836291	,0                   ,0.530073448183213));
    _colorMap.push_back(glm::vec3(0.20914777974203	,0                   ,0.549489950414628));
    _colorMap.push_back(glm::vec3(0.232081027724599	,0                   ,0.565859652831004));
    _colorMap.push_back(glm::vec3(0.254376493411996	,0                   ,0.580276377963287));
    _colorMap.push_back(glm::vec3(0.276334039945701	,0                   ,0.593292780793618));
    _colorMap.push_back(glm::vec3(0.298114528924529	,0                   ,0.605189675391704));
    _colorMap.push_back(glm::vec3(0.319801601946233	,0                   ,0.616096708548274));
    _colorMap.push_back(glm::vec3(0.341447476870458	,0                   ,0.626048229426296));
    _colorMap.push_back(glm::vec3(0.363060872040296	,0                   ,0.635018727717623));
    _colorMap.push_back(glm::vec3(0.384643119406371	,0                   ,0.642940960039945));
    _colorMap.push_back(glm::vec3(0.406174398785712	,0                   ,0.64972044764354 ));
    _colorMap.push_back(glm::vec3(0.427628591987727	,0                   ,0.655246288682272));
    _colorMap.push_back(glm::vec3(0.448972198651563	,0                   ,0.659401089810452));
    _colorMap.push_back(glm::vec3(0.470160954779415	,0                   ,0.662067731238469));
    _colorMap.push_back(glm::vec3(0.491151383136547	,0                   ,0.663139007200098));
    _colorMap.push_back(glm::vec3(0.511893299310341	,0                   ,0.66252659304917 ));
    _colorMap.push_back(glm::vec3(0.532338408070102	,0                   ,0.660171480052806));
    _colorMap.push_back(glm::vec3(0.552438460485469	,0                   ,0.656044166075816));
    _colorMap.push_back(glm::vec3(0.572145603639558	,0                   ,0.650164902255845));
    _colorMap.push_back(glm::vec3(0.591422012109857	,0                   ,0.64259536383122 ));
    _colorMap.push_back(glm::vec3(0.610233627834076	,0                   ,0.633444283949038));
    _colorMap.push_back(glm::vec3(0.628551358039902	,0.031974231567611	,0.622864129615194));
    _colorMap.push_back(glm::vec3(0.646360107273743	,0.065243629915767	,0.611040149614806));
    _colorMap.push_back(glm::vec3(0.663652385843576	,0.092630373137044	,0.598183576690168));
    _colorMap.push_back(glm::vec3(0.680428423102558	,0.117135839713515	,0.584515910928952));
    _colorMap.push_back(glm::vec3(0.696699030479482	,0.139963198621606	,0.570249324366222));
    _colorMap.push_back(glm::vec3(0.712479399760442	,0.16168685918189	,0.555582982556407));
    _colorMap.push_back(glm::vec3(0.727789223040996	,0.182628427996449	,0.540684943535736));
    _colorMap.push_back(glm::vec3(0.742651976433697	,0.202989345213013	,0.525699393529855));
    _colorMap.push_back(glm::vec3(0.757089500307298	,0.222909206377555	,0.510733571800661));
    _colorMap.push_back(glm::vec3(0.771125472669343	,0.242494507086229	,0.495861210821004));
    _colorMap.push_back(glm::vec3(0.784777785478172	,0.261826401691546	,0.48113802685856 ));
    _colorMap.push_back(glm::vec3(0.798063325317957	,0.280977467292973	,0.466581615909416));
    _colorMap.push_back(glm::vec3(0.810991665958445	,0.300005887932574	,0.452201160975367));
    _colorMap.push_back(glm::vec3(0.823570424202093	,0.318969193158374	,0.437983884052625));
    _colorMap.push_back(glm::vec3(0.835797559644867	,0.337915735993297	,0.423911184366801));
    _colorMap.push_back(glm::vec3(0.847668035438685	,0.356895301538528	,0.409951687480209));
    _colorMap.push_back(glm::vec3(0.859168695371703	,0.375951087391978	,0.39607716545887 ));
    _colorMap.push_back(glm::vec3(0.870281255472912	,0.395125934532842	,0.382250454343326));
    _colorMap.push_back(glm::vec3(0.880981911101541	,0.41446173521972	,0.36843969564171 ));
    _colorMap.push_back(glm::vec3(0.89124046719552	,0.433997182425789	,0.354613570201905));
    _colorMap.push_back(glm::vec3(0.901022244901873	,0.453768637360678	,0.340744645646221));
    _colorMap.push_back(glm::vec3(0.910286352711236	,0.473806999379604	,0.326808780698883));
    _colorMap.push_back(glm::vec3(0.918991645887924	,0.494148748417449	,0.312787429930677));
    _colorMap.push_back(glm::vec3(0.927091752826398	,0.514815903601009	,0.298664716578196));
    _colorMap.push_back(glm::vec3(0.934538762199782	,0.535834277320537	,0.2844283398551  ));
    _colorMap.push_back(glm::vec3(0.941282424411384	,0.5572218264946     ,0.270074377537449));
    _colorMap.push_back(glm::vec3(0.947271889081135	,0.57899667584166	,0.255600776071119));
    _colorMap.push_back(glm::vec3(0.95245307357068	,0.601172811062166	,0.241014873127543));
    _colorMap.push_back(glm::vec3(0.956773644962289	,0.623756760430048	,0.226342411414271));
    _colorMap.push_back(glm::vec3(0.960171721130576	,0.646762541870732	,0.211611844494381));
    _colorMap.push_back(glm::vec3(0.962590811342094	,0.670190912889062	,0.196898627390582));
    _colorMap.push_back(glm::vec3(0.963968249129797	,0.694046927317866	,0.182315077604744));
    _colorMap.push_back(glm::vec3(0.964235662782577	,0.718330490798972	,0.168044660817245));
    _colorMap.push_back(glm::vec3(0.963323491108128	,0.743042902305134	,0.154368213794057));
    _colorMap.push_back(glm::vec3(0.96115906854414	,0.768180546101711	,0.14170774443206 ));
    _colorMap.push_back(glm::vec3(0.957656027088606	,0.793742561376316	,0.130635096640637));
    _colorMap.push_back(glm::vec3(0.952749777876738	,0.819713553670359	,0.1219110029042  ));
    _colorMap.push_back(glm::vec3(0.946359929183233	,0.846081977729274	,0.11629705915755 ));
    _colorMap.push_back(glm::vec3(0.938428369965401	,0.872828862778365	,0.114283574624139));
    _colorMap.push_back(glm::vec3(0.92892234497138	,0.899925232813487	,0.115531140128853));
    _colorMap.push_back(glm::vec3(0.917891525763134	,0.927326351154526	,0.117980747432062));
    _colorMap.push_back(glm::vec3(0.905660936999817	,0.954926200744405	,0.115310310744978));
    _colorMap.push_back(glm::vec3(0.894058310302958	,0.98225357930478	,0.081068765570473));

    GLinit();
}

}

#endif
