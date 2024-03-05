#include "particlecollision.h"
collision::collision(double _e_t, double _e_n, double _delta_t,double _gravity)
{
    this->delta_t=_delta_t;
    this->e_n=_e_n;
    this->e_t=_e_t;
    Tc=10*this->delta_t;
    gravity=_gravity*VECTOR(0,0,-1);
    printf("et,en,dt,Tc,g = %.3f %.3f %.3f %.3f %.3f\n",this->e_t,this->e_n,this->delta_t,Tc,_gravity);
}

void collision::setBoundary(BoundaryBox _box)
{
    box=_box;
}

void collision::setstep(int _step)
{
    step=_step;
}

void collision::setnparticles(llu _nparticles)
{
    nparticles=_nparticles;
    MyParticles.resize(_nparticles);
}

void collision::initializeparticles()
{
    for(llu i=0;i<nparticles;i++)
    {
        MyParticles[i].X=VECTOR(0.5,0.5,0.5);
        MyParticles[i].a=VECTOR(0,0,-1);
        MyParticles[i].density=1;
        MyParticles[i].radius=0.1;
        /*
        在这里初始化初始位置，半径，密度
        
                
        */
        double radius = MyParticles[i].radius; 
        MyParticles[i].volume = PAI * 4.0 / 3.0 * radius * radius * radius; 
  
        // 计算质量  
        double density = MyParticles[i].density; 
        MyParticles[i].mass = density * MyParticles[i].volume; 
    }
}

void collision::initializeACTM()
{
    MyACTMs.resize(nparticles);  
    // 对于每个内部的vector，也进行resize  
    for (llu i = 0; i < nparticles; ++i) {  
        MyACTMs[i].resize(nparticles + 1);  
    } 
}

void collision::CauculateCoefACTM()
{ 
    // 多出来那一维度是用来储存和壁面的碰撞系数，
    // 其实把和壁面那一维省略掉用小球和自己的碰撞系数代替也可
    double M_eff;//effective mass
    for(llu i=0;i<nparticles;i++)
    {
        for(llu j=0;j<=nparticles;j++)
        {
            if(i==j) continue;
            if(j==nparticles) //处理壁面
                M_eff=MyParticles[i].mass;
            else
            {
                double m_p,m_q;
                m_p=MyParticles[i].mass;
                m_q=MyParticles[j].mass;
                M_eff=m_p*m_q/(m_p+m_q);
            }
            MyACTMs[i][j].k_n=M_eff*(PAI*PAI+log(e_n)*log(e_n))/(Tc*Tc);
            MyACTMs[i][j].k_t=M_eff*(PAI*PAI+log(e_n)*log(e_n))/(Tc*Tc); //存疑
            MyACTMs[i][j].d_n=-2*M_eff*log(e_n)/Tc;
            MyACTMs[i][j].d_t=-2*M_eff*log(e_t)/Tc;
        }
    }
}

void collision::CalculateForceandMoment()
{
    llu i,j;
    double density,volume,mass,radius;
    VECTOR X,V;
    VECTOR F(0,0,0),M(0,0,0);//转动和平动方程的右边
    double distance;
    VECTOR g;//relative velocity between particle center
    VECTOR g_cp;//relative velocity between contact point
    VECTOR g_ncp,g_tcp;//normal,tangental component
    VECTOR UnitNormalVector;
    for(i=0;i<nparticles;i++)
    {
        F=VECTOR(0,0,0);
        M=VECTOR(0,0,0);                
        density = MyParticles[i].density;
        volume = MyParticles[i].volume;
        mass = MyParticles[i].mass;
        radius = MyParticles[i].radius;        
        X=MyParticles[i].X;
        V=MyParticles[i].V;       
        F=F+mass*gravity; //受到重力
        //F=F-    计算浮力
        /*
        计算IB力，和IB力矩

        */
              
       for(j=0;j<nparticles;j++)//检验i球和j球的距离
       {
        if(i==j) continue;
        distance=abs(X-MyParticles[j].X)-radius-MyParticles[j].radius;
        if(distance<0.0)
        {
            UnitNormalVector=1.0/abs(X-MyParticles[j].X)*(MyParticles[j].X-X);//高危行径 方向问题
            g=V-MyParticles[j].V;
            g_cp=CalculateRelativeV(MyParticles[i],MyParticles[j]);
            g_ncp=(g_cp*UnitNormalVector)*UnitNormalVector;
            g_tcp=g_cp-g_ncp;
            //先算法向力
            F=F+(-MyACTMs[i][j].k_n*distance*UnitNormalVector \
                 -MyACTMs[i][j].d_n*g_ncp);
            //处理切向力，存疑
            F=F+(-MyACTMs[i][j].k_t*distance*UnitNormalVector \
                 -MyACTMs[i][j].d_t*g_tcp); 
            //先不处理力矩

        }
        /*
        润滑力
        else if(distance<    &&distance>0.0)
        */
       }
       //壁面处理
       {
        if (box.left.exists) 
        {  
            distance=abs(box.left.coordinate - X.x)-radius;
            if(distance<0)
            {
                F=F+(-MyACTMs[i][nparticles].k_n*distance*VECTOR(1,0,0) \
                 -MyACTMs[i][nparticles].d_n*(V));
            }  
        }
        if (box.right.exists) 
        {  
            distance=abs(box.right.coordinate - X.x)-radius;
            if(distance<0)
            {
                F=F+(-MyACTMs[i][nparticles].k_n*distance*VECTOR(-1,0,0) \
                 -MyACTMs[i][nparticles].d_n*(V));
            }  
        }
        if (box.bottom.exists) 
        {  
            distance=abs(box.bottom.coordinate - X.z)-radius;
            if(distance<0)
            {
            if(Isco==false)
            {
                Isco=true;
                stepseparate=10+nowstep;
                printf("velocityin\t%.5f\n",MyParticles[i].V.z);
                proportion=MyParticles[i].V.z;
            }
                F=F+(-MyACTMs[i][nparticles].k_n*distance*VECTOR(0,0,1) \
                 -MyACTMs[i][nparticles].d_n*(V));
            }  
        }
        if (box.top.exists) 
        {  
            distance=abs(box.top.coordinate - X.z)-radius;
            if(distance<0)
            {
                F=F+(-MyACTMs[i][nparticles].k_n*distance*VECTOR(0,0,-1) \
                 -MyACTMs[i][nparticles].d_n*(V));
            }  
        }
        if (box.front.exists) 
        {  
            distance=abs(box.front.coordinate - X.y)-radius;
            if(distance<0)
            {
                F=F+(-MyACTMs[i][nparticles].k_n*distance*VECTOR(0,1,0) \
                 -MyACTMs[i][nparticles].d_n*(V));
            }

        }
        if (box.back.exists) 
        {  
            distance=abs(box.left.coordinate - X.y)-radius;  
            if(distance<0)
            {
                F=F+(-MyACTMs[i][nparticles].k_n*distance*VECTOR(0,-1,0) \
                 -MyACTMs[i][nparticles].d_n*(V));
            }
        }
        }
       //修改i球方程的右端项
       MyParticles[i].force=F;
       MyParticles[i].moment=M;
    }
}

void collision::ParticleStatusChange()
{
    for(llu i=0;i<nparticles;i++)
    {
        if(nowstep==stepseparate)
        {
            proportion=MyParticles[i].V.z/proportion;
            Isco=false;
            printf("velocityout\t%.5f\n",MyParticles[i].V.z);
            printf("in/out=%.5f actual\\e_n\n",abs(proportion));
        }
        MyParticles[i].a=MyParticles[i].force*(1.0/MyParticles[i].mass);
        MyParticles[i].V=MyParticles[i].V+MyParticles[i].a*delta_t;
        MyParticles[i].X=MyParticles[i].X+MyParticles[i].V*delta_t;
    }
}

void collision::output()
{
    //先只输出纵坐标
    for(llu i=0;i<nparticles;i++)
    ss<<nowstep*delta_t<<" "<<fixed<<setprecision(4)<<MyParticles[i].X.z<<endl;
}

void BoundaryBox::setBoundary(bool leftExists, double leftCoord,  
                     bool rightExists, double rightCoord,  
                     bool bottomExists, double bottomCoord,  
                     bool topExists, double topCoord,  
                     bool frontExists, double frontCoord,  
                     bool backExists, double backCoord)
{
    left = Boundary(leftExists, leftCoord);  
    right = Boundary(rightExists, rightCoord);  
    bottom = Boundary(bottomExists, bottomCoord);  
    top = Boundary(topExists, topCoord);  
    front = Boundary(frontExists, frontCoord);  
    back = Boundary(backExists, backCoord);    
}

VECTOR CalculateRelativeV(const Particle &a, const Particle &b)
{
   // 计算两个球心之间的单位向量 n  
    VECTOR n = (b.X - a.X).Unit(); // 假设 X 表示位置  
  
    // 计算两个小球之间的相对速度 g  
    VECTOR g = a.V - b.V;  
  
    // 计算 R_a * (omega_a ^ n) 和 R_b * (omega_b ^ n)  
    VECTOR temp1 = a.omega ^ n;  
    VECTOR temp2 = b.omega ^ n;  
    VECTOR result1 = a.radius * temp1;  
    VECTOR result2 = b.radius * temp2;  
  
    // 计算最终结果  
    VECTOR ans = g + result1 + result2;  
  
    return ans;     
}
