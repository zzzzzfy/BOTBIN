#ifndef	CONFIG_H
#define CONFIG_H
#include <string>

using namespace std;
class Config {
public:
	
	string action_ = "";
	string prefix_ = "/home/user_name/dataset";
	string dataset_ = "twitter.txt";
	string dest_ = "/home/user_name/dataset/twitter";
	string src_ = "/home/user_name/dataset/out.twitter";
	string operation_ = "";
	double eps = 0.0;
	int miu = 0;
	int n_;
	void set_file(int _file) {

		// need to set avg & M in calc_k to 
 		if (_file == 0) {
			dest_ = "/home/user_name/dataset/twitter/twitter";
			src_ = "/home/user_name/dataset/twitter/out.twitter";
			n_ = 41652230;
		}
		if (_file == 1) {
			dest_ = "/home/user_name/dataset/pokec/pokec";
			src_ = "/home/user_name/dataset/pokec/out.pokec";
			n_ = 1632803;
		}
		if (_file == 2) {
			dest_ = "/home/user_name/dataset/skitter/skitter";
			src_ = "/home/user_name/dataset/skitter/out.skitter";
			n_ = 1696415;
		}
		if (_file == 3) {
			dest_ = "/home/user_name/dataset/friend/friend";
			src_ = "/home/user_name/dataset/friend/out.friend";
			n_ = 65608367;
		}
		if (_file == 4) {
			dest_ = "/home/user_name/dataset/webbase/webbase";
			src_ = "/home/user_name/dataset/webbase/out.webbase";
			n_ = 115554441;
		}
		if (_file == 5) {
			dest_ = "/home/user_name/dataset/sina/sina";
			src_ = "/home/user_name/dataset/sina/out.sina";
			n_ = 58655849;
		}//sina
		if (_file == 6) {
			dest_ = "/home/user_name/dataset/web12/web12";
			src_ = "/home/user_name/dataset/web12/out.web12";
			n_ = 90320661;
		}
		if (_file == 7) {
			dest_ = "/home/user_name/dataset/brain/brain";
			src_ = "/home/user_name/dataset/brain/out.brain";
			n_ = 784262;
		}
		if (_file == 8) {
			dest_ = "/home/user_name/dataset/orkut/orkut";
			src_ = "/home/user_name/dataset/orkut/out.orkut";
			n_ = 3072441;
		}
		//youtube
		if (_file == 9) {
			dest_ = "/home/user_name/dataset/youtube/youtube";
			src_ = "/home/user_name/dataset/youtube/out.youtube";
			n_ = 3223585;
		}
		if (_file == 10) {
			dest_ = "/home/user_name/dataset/flickr/flickr";
			src_ = "/home/user_name/dataset/flickr/out.flickr";
			//n_ = 41652230;
		}
		if (_file == 11) {
			//pp minier
			dest_ = "/home/user_name/dataset/pp/pp";
			src_ = "/home/user_name/dataset/pp/out.pp";
			n_ = 8254696;
		}
		//lj
		if (_file == 12) {
			dest_ = "/home/user_name/dataset/lj/lj";
			src_ = "/home/user_name/dataset/lj/out.lj";
			n_ = 4847571;
		}
		//topcat
		if (_file == 13) {
			dest_ = "/home/user_name/dataset/topcat/topcat";
			src_ = "/home/user_name/dataset/topcat/out.topcat";
			n_ = 1791489;
		}
	}
};

#endif