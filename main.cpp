/*

Edoardo Carpita
HPC projet 1 - travail d'amelioration sur la librarie PICCANTE v0.5 stable

 Note: l'ancienne version de la librarie originaire Eigen du projet est encore disponible dans les sources sous include/external/Eigen(old)

*/

#include "include/piccante.hpp"

int main(int argc, char *argv[])
{

    printf("Adding file names to the merger...");
    pic::HDRMerger merger;

    for(int i = 0; i < 7; i++) {
        std::string name = "../data/input/stack/stack_room_exp_" + pic::fromNumberToString(i) + ".jpg";
        merger.addFile(name);
    }

    printf("\nOk one\n");

    printf("Merging LDR images into an HDR image...");
    pic::Image *imgOut = merger.execute(NULL);
    printf("\nOk two\n");

    if(imgOut != NULL) {
        if(imgOut->isValid()) {
            imgOut->Write("../data/output/image_debevec_crf.hdr");
            pic::Image *imgTmo = pic::ReinhardTMO::executeGlobal1(imgOut, NULL);
            imgTmo->Write("../data/output/image_debevec_crf_tmo.png", pic::LT_NOR_GAMMA);
            delete imgTmo;
            delete imgOut;
        }
    }
    printf("\nOk done\n");

    return 0;
}
