#pragma once
#include "fecore_api.h"
#include "DumpStream.h"
#include "FEItemList.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
// Forward declarations
class FEMesh;
class FENode;

//-----------------------------------------------------------------------------
//! Defines a node set
//
class FECORE_API FENodeSet : public FEItemList
{
public:
	FENodeSet();
	FENodeSet(FEMesh* pm);
	FENodeSet(const FENodeSet& ns);

	void operator = (const FENodeSet& ns);

	void create(int n);

	void add(int n);

	void add(const std::vector<int>& ns);

	void add(const FENodeSet& ns);

	int size() const { return (int)m_Node.size(); }

	int& operator [] (int i) { return m_Node[i]; }

	const int& operator [] (int i) const { return m_Node[i]; }

	void SetID(int n) { m_nID = n; }
	int GetID() const { return m_nID; }

	std::vector<int>& GetNodeList() { return m_Node; }

	FENode* Node(int i);
	const FENode* Node(int i) const;

	void Serialize(DumpStream& ar);

	FEMesh* GetMesh() { return m_pmesh; }

protected:
	int					m_nID;		//!< ID of nodeset
	std::vector<int>	m_Node;		//!< list of nodes

protected:
	FEMesh*	m_pmesh;	//!< pointer to parent mesh
};