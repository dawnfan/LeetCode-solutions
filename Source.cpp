#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <sstream>
#include <queue>

using namespace std;

#define max(a, b) a<b?b:a

// Definition for a binary tree node.
struct TreeNode {
	int val;
	TreeNode *left;
	TreeNode *right;
	TreeNode(int x) : val(x), left(NULL), right(NULL) {}
};


//Definition for singly-linked list.
struct ListNode {
	int val;
	ListNode *next;
	ListNode(int x) : val(x), next(NULL) {}
};


class Solution {
public:
	int hammingDistance(int x, int y) {
		int xor_val = x ^ y;
		return countBits(xor_val);
	}
	int countBits(int n) {
		// count the 1-bits in a number
		n = (n & 0x55555555) + ((n >> 1) & 0x55555555);
		n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
		n = (n & 0x0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f);
		n = (n & 0x00ff00ff) + ((n >> 8) & 0x00ff00ff);
		n = (n & 0x0000ffff) + ((n >> 16) & 0x0000ffff);
		return n;
	}

	int findComplement(int num) {
		int out_num = 0;
		for (int i = -1; num; num >>= 1)
		{
			if (num & 1 > 0)
			{
				i += num & 1;
				continue;
			}
			i++;
			out_num += 1 << i;
		}
		return out_num;
	}

	bool checkInCode(string w) {
		set<char> codeset[3];
		char code1[] = { 'q', 'w', 'e', 'r', 't', 'y', 'u', 'i', 'o', 'p', 'Q', 'W', 'E', 'R', 'T', 'Y', 'U', 'I', 'O', 'P' };
		codeset[0] = set<char>(code1, code1 + 20);
		char code2[] = { 'a', 's', 'd', 'f', 'g', 'h', 'j', 'k', 'l', 'A', 'S', 'D', 'F', 'G', 'H', 'J', 'K', 'L' };
		codeset[1] = set<char>(code2, code2 + 18);
		char code3[] = { 'z', 'x', 'c', 'v', 'b', 'n', 'm', 'Z', 'X', 'C', 'V', 'B', 'N', 'M' };
		codeset[2] = set<char>(code3, code3 + 14);

		bool isIncode = true;
		for (int i = 0; i < 3; i++)
		{
			isIncode = true;
			set<char> cur_code = codeset[i];
			for (int j = 0; j < w.size(); j++)
			{
				char tt = w[j];
				if (cur_code.find(tt) == cur_code.end())
				{
					isIncode = false;
					break;
				}
			}
			if (isIncode)
			{
				break;
			}
		}
		return isIncode;
	}
	vector<string> findWords(vector<string>& words) {
		vector<string> out_words;
		for (vector<string>::iterator w = words.begin(); w != words.end(); w++)
		{
			if (checkInCode(*w))
			{
				out_words.push_back(*w);
			}
		}
		return out_words;
	}

	vector<string> fizzBuzz(int n) {
		vector<string> outwords;
		for (int i = 1; i < n + 1; i++)
		{
			// check multiples 
			if (i % 3 == 0 && i % 5 == 0)
			{
				outwords.push_back("FizzBuzz");
			}
			else if (i % 3 == 0)
			{
				outwords.push_back("Fizz");
			}
			else if (i % 5 == 0)
			{
				outwords.push_back("Buzz");
			}
			else
			{
				stringstream sstream;
				sstream << i;
				outwords.push_back(sstream.str());
			}
		}
		return outwords;
	}

	string reverseString(string s) {
		int cur_size = s.size();
		string out_s = "";
		for (int i = cur_size - 1; i > -1; i--)
		{
			out_s += s[i];
		}
		return out_s;
	}

	vector<int> nextGreaterElement(vector<int>& findNums, vector<int>& nums) {
		int num_tol = nums.size();
		vector<int> out_nums;
		for (int i = 0; i < findNums.size(); i++)
		{
			int cur_out = -1;
			int cur_pos = num_tol + 1;
			for (int j = 0; j < num_tol; j++)
			{
				if (findNums[i] == nums[j])
				{
					cur_pos = j;
				}
				if (findNums[i] < nums[j] && cur_pos < j)
				{
					cur_out = nums[j];
					break;
				}
			}
			out_nums.push_back(cur_out);
		}
		return out_nums;
	}

	vector<double> averageOfLevels(TreeNode* root) {
		int cur_num = 0;
		int next_num = 0;
		vector<TreeNode*> cur_list;
		vector<double> out_list;
		double cur_sum = 0.0;
		// BFS the tree
		cur_list.push_back(root);
		cur_num = 1;
		while (!cur_list.empty())
		{
			for (int t = 0; t < cur_num; t++)
			{
				TreeNode* cur_node = cur_list[t];
				cur_sum += cur_node->val;
				if (cur_node->left != NULL)
				{
					next_num++;
					cur_list.push_back(cur_node->left);
				}
				if (cur_node->right != NULL)
				{
					next_num++;
					cur_list.push_back(cur_node->right);
				}
			}
			cur_sum /= cur_num;
			out_list.push_back(cur_sum);
			cur_sum = 0.0;
			cur_list.erase(cur_list.begin(), cur_list.begin() + cur_num);
			cur_num = next_num;
			next_num = 0;
		}

		return out_list;
	}

	TreeNode* mergeNodes(TreeNode* n1, TreeNode* n2) {
		TreeNode* out_node;
		int out_val = 0;
		if (n1 == NULL && n2 == NULL) return NULL;
		if (n1 != NULL)
		{
			out_val += n1->val;
		}
		if (n2 != NULL)
		{
			out_val += n2->val;
		}
		out_node = new TreeNode(out_val);

		// recurse on children
		if (n1 != NULL && n2 == NULL)
		{
			out_node->left = mergeNodes(n1->left, NULL);
			out_node->right = mergeNodes(n1->right, NULL);
		}
		else if (n1 == NULL && n2 != NULL)
		{
			out_node->left = mergeNodes(n2->left, NULL);
			out_node->right = mergeNodes(n2->right, NULL);
		}
		else
		{

			out_node->left = mergeNodes(n1->left, n2->left);
			out_node->right = mergeNodes(n1->right, n2->right);
		}

		return out_node;
	}
	TreeNode* mergeTrees(TreeNode* t1, TreeNode* t2) {
		TreeNode* out_tree;
		out_tree = mergeNodes(t1, t2);
		return out_tree;
	}

	string reverseWords(string s) {
		string out_s;
		int leng_num = s.length();
		string temp_s;
		for (int i = 0; i < leng_num + 1; i++)
		{
			if (s[i] != ' ' && i < leng_num)
			{
				temp_s.insert(0, 1, s.at(i));
			}
			else
			{
				if (i < leng_num) temp_s.append(1, ' ');
				out_s.append(temp_s);
				temp_s.clear();
			}
		}
		return out_s;
	}

	TreeNode* trimBST(TreeNode* root, int L, int R) {
		TreeNode* temp_n;
		if (root == NULL) return NULL;
		if (root->val < L)
		{
			temp_n = trimBST(root->right, L, R);
		}
		else if (root->val > R)
		{
			temp_n = trimBST(root->left, L, R);
		}
		else
		{
			root->right = trimBST(root->right, L, R);
			root->left = trimBST(root->left, L, R);
			temp_n = root;
		}
		return temp_n;
	}

	vector<vector<int>> matrixReshape(vector<vector<int>>& nums, int r, int c) {
		int old_r = nums.size();
		int old_c = nums[0].size();

		if (old_r*old_c != r*c) return nums;

		vector<vector<int> > new_nums;
		for (int ri = 0; ri < r; ri++)
		{
			vector<int> cur_nums_r;
			for (int ci = 0; ci < c; ci++)
			{
				int cur_nn = ri*c + ci;
				cur_nums_r.push_back(nums[cur_nn / old_c][cur_nn%old_c]);
			}
			new_nums.push_back(cur_nums_r);
		}

		return new_nums;
	}

	int islandPerimeter(vector<vector<int>>& grid) {
		int r = grid.size();
		int c = grid[0].size();

		int out_per = 0;
		for (int ri = 0; ri < r; ri++)
		{
			for (int ci = 0; ci < c; ci++)
			{
				int cur_n = grid[ri][ci];
				if (cur_n == 1)
				{
					out_per += (ri == 0) + (ri == (r - 1));
					out_per += (ci == 0) + (ci == (c - 1));
				}
				if (ri > 0)
				{
					if (grid[ri - 1][ci] != cur_n) out_per++;
				}
				if (ci > 0)
				{
					if (grid[ri][ci - 1] != cur_n) out_per++;
				}
			}
		}
		return out_per;
	}

	int singleNumber(vector<int>& nums) {
		int out_n = 0;
		int nn = nums.size();
		for (int i = 0; i < nn; i++)
		{
			out_n = out_n ^ nums[i];
		}
		return out_n;
	}

	int maxDepth(TreeNode* root) {
		if (root == NULL) return 0;
		return max(maxDepth(root->left), maxDepth(root->right)) + 1;
	}

	TreeNode* invertTree(TreeNode* root) {
		if (root == NULL)
		{
			return root;
		}
		TreeNode* r_left = root->left;
		TreeNode* r_right = root->right;
		root->left = r_right;
		root->right = r_left;
		invertTree(r_right);
		invertTree(r_left);
		return root;
	}

	bool detectCapitalUse(string word) {
		bool isCorrect = true;
		bool first_isCapital = (word[0] < 65 || word[0] > 90) ? false : true;
		int nn = word.size();
		if (nn > 1)	first_isCapital &= (word[1] < 65 || word[1] > 90) ? false : true;
		if (first_isCapital)
		{
			for (int i = 1; i < nn; i++)
			{
				if (word[i] < 65 || word[i] > 90)	return false;
			}
		}
		else
		{
			for (int i = 1; i < nn; i++)
			{
				if (word[i] >= 65 && word[i] <= 90)	return false;
			}
		}
		return isCorrect;
	}

	int addDigits(int num) {
		// formula - e.g. 38=>3+8=11,11=>1+1=2
		return 1 + (num - 1) % 9;
	}

	vector<int> findDisappearedNumbers(vector<int>& nums) {
		int nn = nums.size();
		vector<int> out_num(nn, 0);
		for (int i = 0; i < nn; i++)
		{
			out_num[nums[i] - 1]++;
		}
		int tt = 0;
		for (int i = 0; i < nn; i++)
		{
			if (out_num.at(tt) != 0)
			{
				out_num.erase(out_num.begin() + tt);
			}
			else
			{
				out_num.at(tt) = i + 1;
				tt++;
			}
		}
		return out_num;
	}
	// more effective 重点放在标记时不改变本来的值
	vector<int> findDisappearedNumbers2(vector<int>& nums) {
		int len = nums.size();
		for (int i = 0; i < len; i++) {
			int m = abs(nums[i]) - 1; // index start from 0
			nums[m] = nums[m] > 0 ? -nums[m] : nums[m];
		}
		vector<int> res;
		for (int i = 0; i < len; i++) {
			if (nums[i] > 0) res.push_back(i + 1);
		}
		return res;
	}

	// bit manipulation to get sum
	int getSum(int a, int b) {
		int out_sum = 0;

		while (b != 0) {
			out_sum = a^b;
			b = (a&b) << 1;
			a = out_sum;
		}
		return out_sum;
	}

	bool checkTofind(vector<int>& tofind, int cur) {
		bool out_c = false;
		if (!tofind.empty())
			for (vector<int>::iterator ti = tofind.begin(); ti != tofind.end(); ti++)
			{
				if (cur == *ti)
				{
					out_c = true;
					break;
				}
			}
		return out_c;
	}
	bool findTarget(TreeNode* root, int k) {
		// 更好的方法：遍历BST得到一个有序的数组，再双向查找
		// 没有考虑BST的特性
		vector<int> tofind;
		queue<TreeNode*> t_list;
		t_list.push(root);
		// include negative situation
		while (!t_list.empty())
		{
			TreeNode* cur_node = t_list.front();
			t_list.pop();
			if (cur_node == NULL) continue;
			int tt_val = cur_node->val;
			t_list.push(cur_node->left);
			t_list.push(cur_node->right);

			if (checkTofind(tofind, tt_val))
			{
				return true;
			}
			tofind.push_back(k - tt_val);
		}
		return false;
	}

	void moveZeroes(vector<int>& nums) {
		int nn = nums.size();
		int toplace = -1;
		int toplace_length = 0;
		int zero_num = 0;
		for (int i = 0; i < nn; i++)
		{
			int cur_num = nums[i];
			if (cur_num == 0)
			{
				toplace = (toplace_length > 0) ? toplace : i;
				toplace_length++;
				zero_num++;
			}
			else
			{
				if (toplace_length > 0)
				{
					nums[toplace] = cur_num;
					toplace++;
				}
			}
		}
		for (int i = nn - zero_num; i < nn; i++)
		{
			nums[i] = 0;
		}
	}

	string tree2str(TreeNode* t) {
		if (t == NULL) return "";

		string s_out = to_string(t->val);
		string s_left = tree2str(t->left);
		string s_right = tree2str(t->right);
		if (s_left.size() > 0 || s_right.size() > 0)
		{
			s_out += "(";
			s_out += s_left;
			s_out += ")";
			if (s_right.size() > 0)
			{
				s_out += "(";
				s_out += s_right;
				s_out += ")";
			}
		}

		return s_out;
	}

	int maxAreaOfIsland(vector<vector<int>>& grid) {
		int count_out = 0;
		int n_row = grid.size();
		int n_col = grid[0].size();
		for (int i = 0; i < n_row; i++)
		{
			for (int j = 0; j < n_col; j++)
			{
				if (grid[i][j] == 1)	count_out = max(count_out, dfsIsland(i, j, grid));
			}
		}
		return count_out;
	}
	int dfsIsland(int i, int j, vector<vector<int>>& grid) {
		if (i >= 0 && i < grid.size() && j >= 0 && j < grid[0].size() && grid[i][j] == 1)
		{
			grid[i][j] = 0;
			return 1 + dfsIsland(i - 1, j, grid) + dfsIsland(i, j - 1, grid) + dfsIsland(i + 1, j, grid) + dfsIsland(i, j + 1, grid);
		}
		return 0;
	}

	// current value only connected with right node.
	// left node should consider right.
	int cur_great = 0;
	void getGreater(TreeNode* root) {
		if (!root) return;
		if (root->right) getGreater(root->right);
		root->val += cur_great;
		cur_great = root->val;
		if (root->left) getGreater(root->left);
	}
	TreeNode* convertBST(TreeNode* root) {
		getGreater(root);
		return root;
	}

	int findShortestSubArray(vector<int>& nums) {
		int nn = nums.size();
		if (nn < 2) return nn;
		map<int, int> map_index, map_count;
		int out_length = nn, out_count = 0;
		for (int i = 0; i < nn; i++)
		{
			if (map_index.find(nums[i]) == map_index.end())
			{
				map_index[nums[i]] = i;
			}
			int cur_count = ++map_count[nums[i]];
			int cur_length = 1 + i - map_index[nums[i]];
			if (cur_count == out_count)	out_length = min(out_length, cur_length);
			else if (cur_count > out_count)
			{
				out_length = cur_length;
				out_count = cur_count;
			}
		}
		return out_length;
	}

	int titleToNumber(string s) {
		int out_num = 0;
		int sn = s.size();
		for (int i = 0; i < sn; i++)
		{
			char cur_s = s[i];
			if (cur_s > 64 && cur_s < 91)
			{
				out_num *= 26;
				out_num += (cur_s - 64);
			}
			else {
				out_num = -1;
				break;
			}
		}
		return out_num;
	}

	int minMoves(vector<int>& nums) {
		// sum up the disparities between the minimum and others.
		int out_move = 0;
		int num_min = INT_MAX;
		for (int i = 0; i < nums.size(); i++)
		{
			num_min = min(num_min, nums[i]);
		}
		for (int i = 0; i < nums.size(); i++)
		{
			out_move += nums[i] - num_min;
		}
		return out_move;
	}

	vector<int> intersection(vector<int>& nums1, vector<int>& nums2) {
		set<int> numset1(nums1.begin(), nums1.end());
		vector<int> out_nums;
		for (int i = 0; i < nums2.size(); i++)
		{
			if (numset1.count(nums2[i]))
			{
				numset1.erase(nums2[i]);
				out_nums.push_back(nums2[i]);
			}
		}
		return out_nums;
	}

	bool canConstruct(string ransomNote, string magazine) {
		bool out_re = true;
		int mn = ransomNote.size();
		for (int i = 0; i < mn; i++)
		{
			int cur_find = magazine.find(ransomNote[i]);
			if (cur_find != -1)
			{
				magazine.erase(magazine.begin() + cur_find);
			}
			else {
				return false;
			}
		}
		return out_re;
	}

	int getLeaf(TreeNode* root, bool isLeft) {
		int out_sum = 0;
		if (root == NULL) return out_sum;
		out_sum += getLeaf(root->right, false);
		out_sum += getLeaf(root->left, true);
		if (root->right == NULL && root->left == NULL && isLeft)
		{
			out_sum += root->val;
		}
		return out_sum;
	}
	int sumOfLeftLeaves(TreeNode* root) {
		int out_sum = 0;
		out_sum = getLeaf(root, false);
		return out_sum;
	}

	bool judgeCircle(string moves) {
		char move_tags[4] = { 'U', 'D', 'L', 'R' };
		int move_counts[4] = { 0 };
		for (int i = 0; i < moves.size(); i++)
		{
			switch (moves[i])
			{
			case 'U':
				move_counts[0]++;
				break;
			case 'D':
				move_counts[1]++;
				break;
			case 'L':
				move_counts[2]++;
				break;
			case 'R':
				move_counts[3]++;
				break;
			default:
				return false;
				break;
			}
		}
		if (move_counts[0] == move_counts[1] && move_counts[2] == move_counts[3]) return true;
		return false;
	}

	// output max/min counting of nums => usage of STD:MAP
	int majorityElement(vector<int>& nums) {
		int out_maj;
		int out_count = 0;

		map<int, int> counts_num;
		for (int i = 0; i < nums.size(); i++)
		{
			if (!counts_num.count(nums[i]))
			{
				counts_num[nums[i]] = 0;
			}
			counts_num[nums[i]]++;
			if (out_maj != nums[i] && counts_num[nums[i]] > out_count)
			{
				out_maj = nums[i];
				out_count = counts_num[nums[i]];
			}
			else if (out_maj == nums[i]) out_count = counts_num[nums[i]];

		}
		return out_maj;
	}

	int findContentChildren(vector<int>& g, vector<int>& s) {
		int out_num = 0;
		// two array sorted
		sort(g.begin(), g.end());
		sort(s.begin(), s.end());
		int i = 0, j = 0;
		while (i < g.size() && j < s.size())
		{
			if (g[i] <= s[j])
			{
				out_num++;
				i++;
			}
			j++;
		}
		return out_num;
	}

	vector<int> twoSum(vector<int>& numbers, int target) {
		int nn = numbers.size();
		vector<int> out_index(2);
		for (out_index[0] = 1, out_index[1] = nn; out_index[0] < out_index[1];)
		{
			if (numbers[out_index[0] - 1] + numbers[out_index[1] - 1] < target)	out_index[0]++;
			else if (numbers[out_index[0] - 1] + numbers[out_index[1] - 1] > target) out_index[1]--;
			else
				break;
		}
		return out_index;
	}

	int firstUniqChar(string s) {
		int out_index = -1;
		map<char, int> s_count;
		for (int i = 0; i < s.size(); i++)
		{
			s_count[s[i]]++;
		}
		for (int i = 0; i < s.size(); i++)
		{
			if (s_count[s[i]] == 1)
			{
				out_index = i;
				break;
			}
		}
		return out_index;
	}

	void getSortedTree(vector<int>& out_arr, TreeNode* root) {
		if (root == NULL) return;
		getSortedTree(out_arr, root->left);
		out_arr.push_back(root->val);
		getSortedTree(out_arr, root->right);
		return;
	}
	int getMinimumDifference(TreeNode* root) {
		int out_min = INT_MAX;
		// sort results
		vector<int> sorted_tree;
		getSortedTree(sorted_tree, root);
		// get related difference
		int nn = sorted_tree.size();
		for (int i = 0; i < nn - 1; i++)
		{
			out_min = min(out_min, sorted_tree[i + 1] - sorted_tree[i]);
		}
		return out_min;
	}

	bool isSameTree(TreeNode* p, TreeNode* q) {
		if (p != NULL && q != NULL) {
			if (p->val == q->val)
			{
				return isSameTree(p->left, q->left) && isSameTree(p->right, q->right);
			}
		}
		else if (p == NULL && q == NULL) { return true; }
		return false;
	}

	void deleteNode(ListNode* node) {
		ListNode* next_n = node->next;
		node->val = next_n->val;
		node->next = next_n->next;
		delete next_n;
	}

	int romanToInt(string s) {
		//Ⅰ（1）、X（10）、C（100）、M（1000）、V（5）、L（50）、D（500）
		unordered_map<char, int> pattern_num;
		pattern_num['I'] = 1;
		pattern_num['X'] = 10;
		pattern_num['C'] = 100;
		pattern_num['M'] = 1000;
		pattern_num['V'] = 5;
		pattern_num['L'] = 50;
		pattern_num['D'] = 500;

		int pre_num = 0;
		int out_num = 0;
		int nn = s.size();
		for (int i = 0; i < nn; i++)
		{
			int cur_num = pattern_num[s[i]];
			out_num += cur_num > pre_num ? -pre_num : pre_num;
			pre_num = cur_num;
		}
		out_num += pre_num;
		return out_num;
	}

	ListNode* reverseNode(ListNode*node, ListNode* pre_node){
		if (node == NULL)	return pre_node;
		ListNode* next_n = node->next;
		node->next = pre_node;
		return reverseNode(next_n, node);
	}
	ListNode* reverseList(ListNode* head) {
		return reverseNode(head, NULL);
	}

	int longestPalindrome(string s) {
		unordered_map<char, int> s_count;
		int nn = s.size();
		for (int i = 0; i < nn; i++)
		{
			s_count[s[i]]++;
		}
		int out_l = 0;
		bool hasSingle = false;
		for (unordered_map<char, int>::iterator s_iter = s_count.begin(); s_iter != s_count.end(); s_iter++)
		{
			if (!hasSingle)
			{
				if (s_iter->second % 2 != 0) hasSingle = true;
			}
			out_l += s_iter->second / 2;
		}
		return hasSingle ? out_l * 2 + 1 : out_l * 2;
	}

	int maximumProduct(vector<int>& nums) {
		sort(nums.begin(), nums.end());
		int nn = nums.size();
		int max_1 = nums[0] * nums[1] * nums[2];
		int max_2 = nums[nn - 1] * nums[nn - 2] * nums[nn - 3];
		return max_1 > max_2 ? max_1 : max_2;
	}

	vector<int> intersect(vector<int>& nums1, vector<int>& nums2) {
		// try sorted vectors. two pointers to match.
		vector<int> *shorter_n = nums1.size() > nums2.size() ? &nums2 : &nums1;
		vector<int> *longer_n = nums1.size() > nums2.size() ? &nums1 : &nums2;
		int sn = shorter_n->size();
		int ln = longer_n->size();
		vector<int> out_inter;
		unordered_map<int, int> count_sn;
		for (int i = 0; i < sn; i++)
		{
			count_sn[shorter_n->at(i)]++;
		}
		for (int i = 0; i < ln; i++)
		{
			int cur_n = longer_n->at(i);
			if (count_sn.count(cur_n) > 0)
			{
				if (count_sn[cur_n] > 0)
				{
					out_inter.push_back(cur_n);
					count_sn[cur_n]--;
				}
			}
		}
		return out_inter;
	}

	int heightAndDiameter(TreeNode* root, int & cur_d){
		if (root == NULL) return 0;
		int left_h = heightAndDiameter(root->left, cur_d);
		int right_h = heightAndDiameter(root->right, cur_d);
		cur_d = max(cur_d, left_h + right_h + 1);
		return max(left_h, right_h) + 1;
	}
	int diameterOfBinaryTree(TreeNode* root) {
		int out_diameter = 0;
		heightAndDiameter(root, out_diameter);
		return out_diameter;
	}

	string convertToBase7(int num) {
		string out_s;
		bool isPositive = num < 0 ? false : true;
		int t = abs(num);
		while (t > 0){
			out_s = to_string(t % 7) + out_s;
			t /= 7;
		}
		return isPositive? out_s : "-"+out_s;
	}

	vector<int> twoSum_noSort(vector<int>& nums, int target) {
		int nn = nums.size();
		vector<int> out_index(2);
		unordered_map<int, int> diff_list;
		for (int i = 0; i < nn; i++)
		{
			if (diff_list.count(target-nums[i]) > 0)
			{
				out_index[0] = diff_list[target - nums[i]];
				out_index[1] = i;
				break;
			}
			diff_list[nums[i]] = i;
		}
		return out_index;
	}
};

void printStringVector(vector<string> strs_in) {
	for (int i = 0; i < strs_in.size(); i++)
	{
		cout << strs_in[i] << endl;
	}
}

int main() {
	Solution leetcode;
	//cout << leetcode.hammingDistance(1, 4) << endl;
	//cout << leetcode.findComplement(5) << endl;

	//string ww[] = { "Hello", "Alaska", "Dad", "Peace" };
	//vector<string> inwords(ww, ww+4);
	//vector<string> outwords = leetcode.findWords(inwords);
	//printStringVector(outwords);

	//printStringVector(leetcode.fizzBuzz(2));

	//cout << leetcode.reverseString("hello") << endl;

	//int nn[] = { 4, 1, 2 };
	//int nnn[] = { 1, 3, 4, 2 };
	//vector<int> fn(nn, nn + 3);
	//vector<int> tn(nnn, nnn + 4);
	//vector<int> out_nums = leetcode.nextGreaterElement(fn, tn);

	//string ts("Let's take LeetCode contest");
	//cout << leetcode.reverseWords(ts) << endl;
	//cout << "before: " << ts << endl;

	//int tt[] = { 1 };
	//vector<int> mylist(tt, tt+1);
	//cout << leetcode.singleNumber(mylist) << endl;

	//int a = 2;
	//int b = 3;
	//cout << leetcode.getSum(a, b) << endl;

	//leetcode.tree2str(NULL);

	//vector<vector<int>> test_t(4);
	//int nn[] = { 1, 1, 0, 0, 0 };
	//test_t[0] = vector<int>(nn, nn + 5);
	//test_t[1] = vector<int>(nn, nn + 5);
	//int nnn[] = { 0, 0, 0, 1, 1 };
	//test_t[2] = vector<int>(nnn, nnn + 5);
	//test_t[3] = vector<int>(nnn, nnn + 5);
	//cout << leetcode.maxAreaOfIsland(test_t) << endl;


	//int nnn[] = { 1, 2, 2, 3, 1 };
	//vector<int> tt(nnn, nnn + 5);
	//cout << leetcode.findShortestSubArray(tt) << endl;

	//cout << leetcode.canConstruct("aa", "aab") << endl;

	int nnn[] = { 3, 2, 1, 4, 5 };
	vector<int> tt(nnn, nnn + 5);
	cout << leetcode.maximumProduct(tt) << endl;

	return 0;
}