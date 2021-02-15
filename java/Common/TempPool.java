
package Common;

public class TempPool {

    int icur=0;
    Object [] objs;
    
    public TempPool(int n, Class c){
        objs=new Object[n];
        try{
            for(int i=0; i<n; i++){
                //objs[i]=c.getConstructor().newInstance(objs);
                objs[i]=c.newInstance();
            }
        }catch(Exception e){
            System.out.println( e );
            objs[-1]=null;
            //System.out.println( "Exeption: TempPool.payback() got different instance "+o+" expected"+objs[i] );
        }
    }
    
    public final Object borrow(){
        Object out = objs[icur];
        icur++;
        return out;
    }
    public final boolean payback(Object o){
        int i=icur-1;
        if( objs[i]!=o ){
            System.out.println( "Exeption: TempPool.payback() got different instance "+o+" expected"+objs[i] );
            objs[-1]=null;
            return false;
        }
        icur=i;
        return true;
    }
    
}
