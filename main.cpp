#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include<algorithm>
#include "Eigen/Core"
#include "Eigen/LU"

#define numOfNodeInElm 3
#define weight 1

using namespace Eigen;
using namespace std;

typedef Matrix<double,2,2> Matrix2_2d;
typedef Matrix<double,3,2> Matrix3_2d;
typedef Matrix<double,2,3> Matrix2_3d;
typedef Matrix<double,3,3> Matrix3_3d;


void export_vtu(const std::string &file, vector<vector<double>> node, vector<vector<int>> element, vector<double> C)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
  fprintf(fp, "</Points>\n");
  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < element[i].size(); j++) fprintf(fp, "%d ", element[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += element[i].size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++) fprintf(fp, "%d\n", 5);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");
  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size();
  fprintf(fp, "</PointData>\n");
  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);
  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[node.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < node.size(); ic++){
    data_d[num] = node[ic][0];
    num++;
    data_d[num] = node[ic][1];
    num++;
    data_d[num] = 0.0;
    num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
  num=0;
  for (int ic = 0; ic < node.size(); ic++){
      data_d[num]   = C[ic];
      num++;
  }
  size=sizeof(double)*node.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
  delete data_d;
  ofs.close();
  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}


int main()
{
    string str;
    ifstream ifs("node.dat");
    //vector<double> t(2) は double t[2]と同じ
    //vector<vector<double>> t(2, vector<double>(2))は double t[2][2]と同じ
    vector<vector<double>> x;
    while(getline(ifs,str)){
        istringstream ss(str); 
        string tmp;
        vector<double> tmp_x;
        for(int j=0; j<3; j++){
            getline(ss, tmp, ' ');
            tmp_x.push_back(stod(tmp));
        }
        x.push_back(tmp_x);
    }
    ifs.close();
    ifs.open("element.dat");
    vector<vector<int>> element;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        vector<int> tmp_element;
        for(int j=0; j<4; j++){
            getline(ss, tmp, ' ');
            if(j==0) continue;
            tmp_element.push_back(stoi(tmp));
        }
        element.push_back(tmp_element);
    }
    ifs.close();


    vector<double> C1(x.size(),0.0);
    double minimum = 100000.0;
    for(int i=0; i<x.size(); i++){
        minimum=min(minimum,x[i][0]);
    }

    ofstream ofs1("boundary_left.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][0]-minimum)<0.000001){

            ofs1 << i << endl;
            C1[i] = 1.0;
        }
    }
    ofs1.close();

    
    vector<double> C2(x.size(),0.0);
    double maximum = 0.01;
    for(int i=0; i<x.size(); i++){
        maximum=max(maximum,x[i][0]);
    }

    ofstream ofs2("boundary_right.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][0]-maximum)<0.000001){

            ofs2 << i << endl;
            C2[i] = 1.0;
        }
    }

    ofs2.close();


    export_vtu("test.vtu", x, element, C1);
    

    ifs.open("boundary_left.dat");
    vector<int> boundary_left;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
        boundary_left.push_back(stoi(tmp));
    }
    ifs.close();

    ifs.open("boundary_right.dat");
    vector<int> boundary_right;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
        boundary_right.push_back(stoi(tmp));
    }
    ifs.close();


    Matrix3_2d dNds;
    Matrix3_2d dNdx;
    Matrix2_2d dsdx;
    Matrix2_2d dxds;
    Matrix2_3d x_current;
    Matrix3_3d Ke;
    MatrixXd K(1301,1301);
    VectorXd U(1301);
    VectorXd R(1301);
    U = VectorXd::Zero(1301);
    R = VectorXd::Zero(1301);
    

    double detJ;

    dNds << -1, -1,
             1, 0,
             0, 1;
    
    for(int i=0; i<x.size(); i++){
      for(int j=0; j<x.size(); j++){
        K(i,j) = 0;
      }
    }

    for(int ic=0; ic<element.size(); ic++){

      for(int i=0; i<numOfNodeInElm; i++){
        for(int j=0; j<2; j++){
          x_current(j,i) = 0;
          x_current(j,i) = x[element[ic][i]][j];
        }
      }

      //initialize
      for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
          dxds(i,j) = 0;
          dsdx(i,j) = 0;
        }
      }

      for(int p=0; p<2; p++){
        for(int i=0; i<2; i++){
          dxds(p,i) = 0;
          for(int j=0; j<numOfNodeInElm; j++){
            dxds(p,i) += dNds(j,i)*x_current(p,j);
          }
        }
      }     

      detJ = dxds.determinant(); 

      dsdx = dxds.inverse();

      for(int p=0; p<3; p++){
        for(int i=0; i<2; i++){
          dNdx(p,i) = 0;
          for(int j=0; j<2; j++){
            dNdx(p,i) += dNds(p,j)*dsdx(j,i);
          }
        }
      }

      for(int p=0; p<numOfNodeInElm; p++){
        for(int q=0; q<numOfNodeInElm; q++){
          Ke(p,q) = 0;
          for(int i=0; i<2; i++){
            Ke(p,q) += dNdx(p,i)*dNdx(q,i)*detJ*weight*weight;
          }
        }
      }
  
      for(int p=0; p<numOfNodeInElm; p++){
        for(int q=0; q<numOfNodeInElm; q++){
          K(element[ic][p],element[ic][q]) += Ke(p,q);
        }
      }
      /*ofstream outputfile("K.txt");
      outputfile<< K ;
      outputfile.close();
      exit(1);*/
    }

    /*ofstream outputfile("K.txt");
    outputfile<< K ;
    outputfile.close();*/
    
    //neuman boundary condition
    for(int i=0; i<boundary_left.size(); i++){
      for(int j=0; j<x.size(); j++){
        K(boundary_left[i],j) = 0;
      }
      K(boundary_left[i],boundary_left[i]) = 1;
      R(boundary_left[i]) = 100;
    }

    for(int i=0; i<boundary_right.size(); i++){
      for(int j=0; j<x.size(); j++){
        K(boundary_right[i],j) = 0;
      }
      K(boundary_right[i],boundary_right[i]) = 1;
      R(boundary_right[i]) = 2;
    }

      FullPivLU<MatrixXd> lu(K);
      U = lu.solve(R);
      ofstream outputfile("U.txt");
      outputfile << U;
      outputfile.close();

      vector<double> p(x.size());
      for(int i=0; i<x.size(); i++){
        p[i]=U(i);
      }
    
      export_vtu("result.vtu",x,element,p);
      
  return 0;
}

