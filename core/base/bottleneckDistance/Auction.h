/// \ingroup base
/// \class ttk::Auction
/// \author Joseph Budin <joseph.budin@polytechnique.edu>

#ifndef _AUCTION_H
#define _AUCTION_H

#ifndef matchingTuple
#define matchingTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::idVertex, dataType>
#endif

#ifndef diagramTuple
#define diagramTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::NodeType, ttk::ftm::idVertex, \
  ttk::ftm::NodeType, dataType, ttk::ftm::idVertex, \
  dataType, float, float, float, dataType, float, float, float>
#endif


#include <cmath>
#include <limits>
#include <iostream>
#include <Debug.h>
#include <PersistenceDiagram.h>
#include <AuctionActor.h>

namespace ttk {
    
  template<typename dataType>
  class Auction : public Debug
  {

    public:

		Auction(int wasserstein) {
            n_bidders_ = 0;
            n_goods_ = 0;
			epsilon_ = 1;
			wasserstein_ = wasserstein;
			delta_lim_ = 0.01;
        };
		~Auction() {};

		dataType run(std::vector<matchingTuple> *matchings);

		
		void BuildAuctionDiagrams(std::vector<diagramTuple> diagram1, std::vector<diagramTuple> diagram2){
			n_bidders_ = diagram1.size();
			n_goods_ = diagram2.size();
			this->setBidders(diagram1);
			this->setGoods(diagram2);
			for(int i=0; i < n_bidders_; i++){
				//Add diagonal goods
				Bidder<dataType>& b = bidders_.get(i);
				Good<dataType> g = Good<dataType>(b.x_, b.y_, true, -b.id_-1);
				g.projectOnDiagonal();
				diagonal_goods_.addGood(g);
			}
			for(int i=0; i < n_goods_; i++){
				//Add diagonal bidders
				Good<dataType>& g = goods_.get(i);
				Bidder<dataType> b = Bidder<dataType>(g.x_, g.y_, true, -g.id_-1);
				b.projectOnDiagonal();
				b.setPositionInAuction(bidders_.size());
				bidders_.addBidder(b);
			}
		}
		
		
		void setBidders(std::vector<diagramTuple> diagram1){
			int d1Size = (int) diagram1.size();
			
			for(int i=0; i<d1Size; i++){
				//Add bidder to bidders
				Bidder<dataType> b = Bidder<dataType>(diagram1[i], i);
				b.setPositionInAuction(bidders_.size());
				bidders_.addBidder(b);
			}
			n_bidders_ = bidders_.size();
		}

		void setGoods(std::vector<diagramTuple> diagram2){
			int d2Size = (int) diagram2.size();
			
			for(int i=0; i<d2Size; i++){
				//Add bidder to bidders
				Good<dataType> g = Good<dataType>(diagram2[i], i);
				goods_.addGood(g);
			}	
			n_goods_ = goods_.size();
		}		
		
		void setEpsilon(dataType epsilon){
			epsilon_ = epsilon;
		}
		
		void initializeEpsilon(){
			dataType max_persistence = 0;
			for(int i=0; i<bidders_.size(); i++){
				Bidder<dataType>& b = bidders_.get(i);
				dataType persistence = b.getPersistence();
				if(persistence>max_persistence){
					max_persistence = persistence;
				}
			}
			
			for(int i=0; i<goods_.size(); i++){
				Good<dataType>& g = goods_.get(i);
				dataType persistence = g.getPersistence();
				if(persistence>max_persistence){
					max_persistence = persistence;
				}
			}
			this.epsilon = 5/4 * pow(max_persistence, wasserstein_);
		}
		
		
		void buildUnassignedBidders(){
			for(int i=0; i<bidders_.size(); i++){
				Bidder<dataType>& b = bidders_.get(i);
				b.setProperty(NULL);
				unassignedBidders_.push_back(i);
			}
		}
		
		
		void reinitializeGoods(){
			for(int i=0; i<goods_.size(); i++){
				Good<dataType>& g = goods_.get(i);
				g.setOwner(-1);
			}
			for(int i=0; i<diagonal_goods_.size(); i++){
				Good<dataType>& g = diagonal_goods_.get(i);
				g.setOwner(-1);
			}
		}
		
		
		dataType getMatchingDistance(){
			dataType d = 0;
			for(int i; i<bidders_.size(); i++){
				Bidder<dataType>& b = bidders_.get(i);
				d += b.cost(b.getProperty(), wasserstein_); 
			}
			return d;
		}
		
		
		dataType getRelativePrecision(){
			dataType d = this->getMatchingDistance();
			if(d<1e-12){
				return 0;
			}
			dataType denominator = d - bidders_.size()*epsilon_;
			if(denominator<=0){
				return 1;
			}
			else{
				return pow(d/denominator, 1/((float)wasserstein_)) - 1;
			}
		}
		
		
		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}
	
    protected:
		int wasserstein_;   // Power in Wassertsein distance (by default set to 2)
		BidderDiagram<dataType>  bidders_;
		GoodDiagram<dataType>  goods_;
		GoodDiagram<dataType> diagonal_goods_;
		std::list<int> unassignedBidders_;
		
		int n_bidders_;
		int n_goods_;
		
		dataType epsilon_;
		dataType delta_lim_;
  };
}

#include <AuctionImpl.h>

#endif

