#ifndef MODULE_H_DEF
#define MODULE_H_DEF

struct Module
{
    std::string name;
    long long minArea;
    long long width, height;
    long long cordX, cordY;     // module左下角座標
    long long rightCordX, rightCordY;   // module右上角座標
    std::vector<std::pair<long long, long long>> vec_allpossibleShape;     // width, height
    double aspectRatio;
    long long minWidth, maxWidth;
    bool isFixed;
};

#endif