void printData()
{
    gather();

    printf("\nData in channel R: \n");
    for(int i = 0; i < mesh->N(); i++)
    {
        printf("%4.2f ",h_data[i].x);
        if(i%mesh->n<0>() == mesh->n<0>()-1) {printf("\n");}
    }

    printf("\nData in channel G: \n");
    for(int i = 0; i < mesh->N(); i++)
    {
        printf("%4.2f ",h_data[i].y);
        if(i%mesh->n<0>() == mesh->n<0>()-1) {printf("\n");}
    }

    printf("\nData in channel B: \n");
    for(int i = 0; i < mesh->N(); i++)
    {
        printf("%4.2f ",h_data[i].z);
        if(i%mesh->n<0>() == mesh->n<0>()-1) {printf("\n");}
    }

    printf("\nData in channel A: \n");
    for(int i = 0; i < mesh->N(); i++)
    {
        printf("%4.2f ",h_data[i].w);
        if(i%mesh->n<0>() == mesh->n<0>()-1) {printf("\n");}
    }
}
template<unsigned int Y>
void printDataChannel()
{
    static_assert((Y < 4), "channel index must be 0, 1, 2 or 3!");

    gather();

    printf("\nData in channel %d: \n",Y);
    for(int i = 0; i < mesh->N(); i++)
    {
        printf("%4.2f ",h_data[i][Y]);
        if(i%mesh->n<0>() == mesh->n<0>()-1) {printf("\n");}
    }
}
