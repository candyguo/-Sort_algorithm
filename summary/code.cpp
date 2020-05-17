#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <map>
#include <set>
#include <queue>
#include <tuple>
#include <chrono>
#include <thread>
#include <mutex>

// 每个数字都应该在数字所对应的位置上面
int find_repeated_num(std::vector<int>& nums) {
    std::vector<int> postions(nums.size(),-1);
    for(int i = 0; i < nums.size(); i++) {
        int num = nums.at(i);
        if(postions[num] != -1) {
            return num;
        } else {
            postions[num] = num;
        }
    }
    return -1;
}

bool is_valid_pair(const char a, const char b) {
    if(a == '(')
        return b == ')';
    if(a == '[')
        return b == ']';
    if(a == '{')
        return b == '}';
    return false;            
}

bool is_valid_bracket(const std::string& str) {
    if(str.empty())
        return true;    
    std::stack<char> bracket_stack;
    bracket_stack.push(str.front());
    for(size_t i = 1; i < str.size(); i++) {
        if(bracket_stack.empty()) {
            bracket_stack.push(str[i]);
            continue;
        }
        auto c = bracket_stack.top();
        if(is_valid_pair(c,str[i])) {
            bracket_stack.pop();
        } else {
            bracket_stack.push(str[i]);
        }
    }
    return bracket_stack.empty();
}

class min_stack {
public:
    min_stack() {}

    void push(int x) {
        s.push(x);
        v.push_back(x);
        order_v = v;
        std::sort(order_v.begin(),order_v.end());
    }

    void pop() {
        s.pop();
        v.pop_back();
        order_v = v;
        std::sort(order_v.begin(),order_v.end());
    }

    int top() {
        return s.top();
    }
    
    int get_min() {
        return order_v.front();
    }

private:
    std::stack<int> s;
    std::vector<int> v;
    std::vector<int> order_v;
};

//也可以维护一个每次推入元素的最小元素栈

class min_stack_2 {
public:
    min_stack_2() {}

    void push(int x) {
        s.push(x);
        if(min_s.empty())
            min_s.push(x);
        else {
            int t = min_s.top();
            min_s.push(std::min(x,t));
        }    
    }

    void pop() {
        s.pop();
        min_s.pop();
    }

    int top() {
        return s.top();
    }
    
    int get_min() {
        return min_s.top();
    }

    std::stack<int> s;
    std::stack<int> min_s;
};

std::vector<int> next_greater_num(const std::vector<int>& nums1,
                                  const std::vector<int>& nums2) {
    std::vector<int> results(nums1.size(),-1);
    for(size_t i = 0; i < nums1.size(); i++) {
        int num = nums1.at(i);
        auto index = std::find(nums2.begin(),nums2.end(),num);
        int start_index = std::distance(nums2.begin(),index);
        for(int j = start_index; j < int(nums2.size()); j++) {
            if(nums2.at(j) > num) {
                results.at(i) = nums2.at(j);
                break; 
            }
        }
    }
    return results;
}

std::string erase_outer_bracket(const std::string& str) {
    // step1: 先做原语化分解
    std::vector<std::string> primative_strings;
    int start_index = 0;
    int len = 0;
    std::stack<char> bracket_stack;
    for(size_t i = 0; i < str.size(); i++) {
        if(bracket_stack.empty() && len != 0) {
            std::string primative_string = str.substr(start_index,len);
            primative_strings.push_back(primative_string);
            start_index = i;
            len = 0;
        }
        if(bracket_stack.empty()) {
            bracket_stack.push(str[i]);
            len++;
            continue;
        }
            
        auto top = bracket_stack.top();
        if(is_valid_pair(top,str[i])) {
            bracket_stack.pop(); 
        } else {
            bracket_stack.push(str[i]);
        }
        len++;
    }
    if(bracket_stack.empty() && len != 0) {
        std::string primative_string = str.substr(start_index,len);
        primative_strings.push_back(primative_string);
    }
    // step2：对分解得到的序列字符串去掉外层括号做合并
    std::string result = "";
    for(auto& e : primative_strings) {
        result += e.substr(1,int(e.size()) - 2);
    }
    return result;
}

int baseball_score(const std::vector<std::string>& ops) {
    std::vector<int> scores;
    for(int i = 0; i < ops.size(); i++) {
        if(ops[i] == "+") {
            int last_one = scores[scores.size() - 1];
            int last_two = scores[scores.size() - 2];
            scores.push_back(last_one + last_two);
        } else if(ops[i] == "D") {
            int last_one = scores[scores.size() - 1];
            scores.push_back(last_one * 2);
        } else if(ops[i] == "C") {
            scores.pop_back();
        } else {
            scores.push_back(atoi(ops[i].c_str()));
        }
    }
    return std::accumulate(scores.begin(),scores.end(),0);
}

std::string remove_repeat(const std::string& s) {
    std::stack<char> unique_stack;
    for(size_t i = 0; i < s.size(); i++) {
        if(unique_stack.empty()) {
            unique_stack.push(s[i]);
            continue;
        }
        auto c = unique_stack.top();
        if(c == s[i])
            unique_stack.pop();
        else
        {
            unique_stack.push(s[i]);
        }
            
    }
    // 重构result
    std::string result = "";
    while(!unique_stack.empty()) {
        result = unique_stack.top() + result;
        unique_stack.pop();
    }
    return result;
}

bool backspace_compare(const std::string& s, const std::string& t) {
    std::stack<char> stack_s;
    std::stack<char> stack_t;
    for(size_t i = 0; i < s.size(); i++) {
        if(s[i] == '#') {
            if(!stack_s.empty()) {
                stack_s.pop();
            }
        } else {
            stack_s.push(s[i]);
        }
    }

    for(size_t i = 0; i < t.size(); i++) {
        if(t[i] == '#') {
            if(!stack_t.empty()) {
                stack_t.pop();
            }
        } else {
            stack_t.push(s[i]);
        }
    }

    while(!stack_s.empty()) {
        if(stack_s.top() != stack_t.top())
            return false;
        stack_s.pop();
        stack_t.pop();    
    }
    return true;
}

int remove_value(std::vector<int>& arr, const int value) {
    //找到相应的元素，原地覆盖
    int index = 0;
    for(auto& e : arr) {
        if(e != value) {
            arr[index++] = e;
        }
    }
    return index;
}

// use hash table to solve it
std::vector<int> twoSum(std::vector<int>& nums, int target) {
    //hash表储存的是数字及其对应的(位置索引+1)
    std::unordered_map<int,int> hash_table;
    for(size_t i = 0; i < nums.size(); i++) {
        int other_number = target - nums[i];
        if(hash_table[other_number] > 0) {
            return {hash_table[other_number] - 1, static_cast<int>(i)};
        } else {
            // 把当前数先加入哈希表
            hash_table[nums[i]] = i + 1;
        }
    }
    std::cout<<"no pair found"<<std::endl;
    return {};
}

int maxEqualRowsAfterFlips(std::vector<std::vector<int>>& matrix) {
    //翻转后的最多行 = 某一行出现的个数 + 其补行出现的个数
    std::map<std::string,int> has; //遍历使用map来进行遍历
    for(size_t i = 0; i < matrix.size(); i++) {
        std::string tmp = "";
        for(size_t j = 0; j < matrix[i].size(); j++) {
            tmp += std::to_string(matrix[i][j]);
        }
        has[tmp] ++;
    }
    int result = 0;
    for(auto iter = has.begin(); iter != has.end(); iter++) {
        std::string tmp = iter->first;
        std::string rev = "";
        for(size_t i = 0; i < tmp.size();i++) {
            if(tmp[i] == '1')
                rev += '0';
            else
                rev += '1';    
        }
        //原本的每一行都对应一个最多行的个数
        result = std::max(result, iter->second + has[rev]);
    }
    return result;        
}

int singleNumber(std::vector<int>& nums) {
    std::map<int,int> has;
    for(size_t i = 0; i < nums.size(); i++) {
        has[nums[i]]++;
    }
    for(auto iter = has.begin(); iter != has.end(); iter++) {
        if(iter->second == 1)
            return iter->first;
    }
    return nums.front();
}

std::vector<int> intersection(std::vector<int>& nums1, std::vector<int>& nums2) {
    //使用哈希来判断数组交集,准备一半(O(N))
    std::unordered_map<int,int> has1;
    for(auto& num : nums1) {
        has1[num]++;
    }
    std::set<int> results_set;
    for(auto& num : nums2) {
        if(has1[num] > 0) {
            results_set.insert(num);
        }
    }
    std::vector<int> results;
    for(auto& num : results_set) {
        results.push_back(num);
    }
    return results;
}

int fourSumCount(std::vector<int>& A, std::vector<int>& B, 
                 std::vector<int>& C, std::vector<int>& D) {
    //准备一半 ，哈希AB和CD的和
    std::map<int,std::vector<std::pair<int,int>>> hash_AB;
    std::map<int,std::vector<std::pair<int,int>>> hash_CD;
    for(size_t i = 0; i < A.size(); i++) {
        for(size_t j = 0; j < B.size(); j++) {
            int sum = A[i] + B[j];
            hash_AB[sum].push_back({i,j});
        }
    }
    for(size_t i = 0; i < C.size(); i++) {
        for(size_t j = 0; j < D.size(); j++) {
            int sum = C[i] + D[j];
            hash_CD[sum].push_back({i,j});
        }
    }
    int result = 0;
    for(auto iter = hash_AB.begin(); iter != hash_AB.end(); iter++) {
        int value = -(iter->first);
        if(hash_CD.count(value) != 0)
            result += iter->second.size() * hash_CD[value].size();
    }
    return result;    
}

class TinyUrl {
public:
    std::string encode(std::string long_url) {
        return long_url;
    }

    std::string decode(std::string short_url) {
        return short_url;
    }
};


int subarraySum(std::vector<int>& nums, int k) {
    // 计算累积分布
    std::map<int,int> accumulate_map;
    accumulate_map[0] = nums.front();
    for(size_t i = 1; i < nums.size(); i++) {
        accumulate_map[i] = nums[i] + accumulate_map[i-1];
    }
    std::map<int,std::vector<int>> accumulate_value;
    for(auto iter = accumulate_map.begin(); iter != accumulate_map.end(); iter++) {
        int index = iter->first;
        int sum = iter->second;
        accumulate_value[sum].push_back(iter->first);
    }
    int result = 0;
    if(accumulate_value.count(k) != 0)
        result += accumulate_value[k].size();
    int start_index = 1;
    for( ; start_index < int(nums.size()); start_index++) {
        int accumulate_sum = k + accumulate_map[start_index - 1];
        if(accumulate_value.count(accumulate_sum) != 0) {
            auto v = accumulate_value[accumulate_sum];
            result += std::count_if(v.begin(),v.end(),[&](int i) {
                return i >= start_index;
            });
        }
            
    }
    return result;
}


// todo use hash
int subarraysDivByK(std::vector<int>& A, int K) {
    std::map<int,int> accumulate_map;
    accumulate_map[0] = A.front();
    for(size_t i = 1; i < A.size(); i++) {
        accumulate_map[i] = A[i] + accumulate_map[i-1];
    }
    std::map<int,std::vector<int>> accumulate_value;
    for(auto iter = accumulate_map.begin(); iter != accumulate_map.end(); iter++) {
        int index = iter->first;
        int sum = iter->second;
        accumulate_value[sum % K].push_back(iter->first);
    }
    int result = 0;
    if(accumulate_value.count(0) != 0)
        result += accumulate_value[0].size();  
    int start_index = 1;
    for(; start_index < A.size(); start_index++) {
        int prenum_sum_mod = accumulate_map[start_index - 1];
        for(int j = start_index; j < A.size();j ++) {
            int value = accumulate_map[j] - prenum_sum_mod;
            if(value % K == 0)
                result++;
        }
    }
    return result;
}

struct TreeNode {
    int val;
    TreeNode* left_child;
    TreeNode* right_child;
    TreeNode(int v) : val(v),left_child(nullptr),right_child(nullptr) {}
    TreeNode(int v, TreeNode* left, TreeNode* right) :
             val(v),left_child(left),right_child(right) {};
};

//递归结构恢复二叉树，哈希表进行元素查询
class FindElements {
public:
    FindElements(TreeNode* root) {
        recover(root,0);
    }
    
    bool find(int target) {
        return has[target] > 0;
    }

    void recover(TreeNode* root,int val) {
        if(root == nullptr)
            return;
        root->val = val;
        has[val]++;
        recover(root->left_child,2*val + 1);
        recover(root->right_child,2*val + 2);
    }

    // default value = 0;
    std::map<int,int> has;
};

void moveZeroes(std::vector<int>& nums) {
    int pos = 0;
    for(auto& num : nums) {
        if(num != 0) {
            nums[pos] = num;
            pos++;
        }
    }
    for(int i = pos; i < nums.size(); i++) {
        nums[i] = 0;
    }    
}

int removeDuplicates(std::vector<int>& nums) {
    auto iter = std::unique(nums.begin(),nums.end());
    return std::distance(nums.begin(),iter);
}

int removeDuplicates_2(std::vector<int>& nums) {
    if(nums.empty())
        return 0;
    int pos = 1;
    int repeat_times = 1;
    for(size_t i = 1; i < nums.size(); i++) {
        if(nums[i] != nums[i-1]) {
            nums[pos++] = nums[i];
            repeat_times = 1;
        } else {
            if(repeat_times == 1) {
                nums[pos++] = nums[i];
                repeat_times = 2;
            } else {
                continue;
            }
        }
    }
    return pos;
}

void reverseString(std::vector<char>& s) {
    std::string a;
    std::reverse(s.begin(),s.end());
}

//翻转索引 i 和 j 之间的数据，包括索引 i 和 j
void reverseStr(std::string& s, int i, int j) {
    int ii = i;
    int jj = j;
    while(ii < jj) {
        std::swap(s[ii],s[jj]);
        ii++;
        jj--;
    }
}

//索引的递增
std::string reverseStr(std::string s, int k) {
    std::string result = s;
    if(k >= result.size())
        reverseStr(result,0,result.size()-1);
    else
        reverseStr(result,0,k-1);
    std::cout<<result<<std::endl;
    int start_index = 2 * k;
    int end_index = start_index + 2 * k - 1;
    while(end_index < result.size()) {
        reverseStr(result,start_index, start_index + k - 1);
        start_index += 2 * k;
        end_index += 2 * k;
    }
    if(start_index + k  >= result.size())
        reverseStr(result,start_index,result.size() - 1);
    else
        reverseStr(result,start_index,start_index + k - 1);
    return result;
}

//无符号整数中1的个数
int hammingWeight(uint32_t n) {
    int cnt = 0;
    while(n > 0) {
        n = n & (n-1);
        cnt++;
    }
    return cnt;   
}

bool isPowerOfTwo(int n) {
    int cnt = 0;
    while(n > 0) {
        cnt += n & 1;
        n = n >> 1;
    }
    return cnt == 1;        
}

std::vector<std::vector<int>> threeSum(std::vector<int>& nums) {
    std::vector<std::vector<int>> results;
    if(nums.size() < 3)
        return results;
    std::map<int,std::vector<std::pair<int,int>>> two_sum_map;
    for(size_t i = 0; i < nums.size(); i++) {
        for(size_t j = i + 1; j < nums.size(); j++) {
            int sum = nums[i] + nums[j];
            two_sum_map[sum].push_back({i,j});
        }
    }
    std::set<std::vector<int>> tmp_results;
    for(size_t i = 0; i < nums.size() - 2; i++) {
        if(two_sum_map.count(-nums[i]) != 0) {
            auto& pairs = two_sum_map[-nums[i]];
            for(auto& pair : pairs) {
                if(pair.first <= i || pair.second <= i)
                    continue;
                std::vector<int> tmp{nums[i],nums[pair.first],nums[pair.second]};
                std::sort(tmp.begin(),tmp.end());
                tmp_results.insert(tmp);   
            }
        }
    }
    for(auto& result : tmp_results) {
        results.push_back(result);
    }  
    return results;
}

std::vector<int> dailyTemperatures(std::vector<int>& T) {
    // std::vector<int> results(T.size(), 0);
    // for(size_t i = 0; i < T.size()-1; i++) {
    //     for(int j = i + 1; j < T.size(); j++) {
    //         if(T[j] > T[i]) {
    //             results[i] = j - i;
    //             break;
    //         }
    //     }
    // }
    // return results;
    
    // 维护一个单调递减栈，有大的数就可以更新结果
    // 递减栈
    std::vector<int> results(T.size(), 0);
    std::stack<std::pair<int,int>> temper_stack; // index and temperatures
    temper_stack.push({0, T.front()});
    for(size_t i = 1; i < T.size(); i++) {
        if(T[i] <= temper_stack.top().second) {
            temper_stack.push({i,T[i]});
            continue;
        }
        while(!temper_stack.empty() && T[i] > temper_stack.top().second) {
            results[temper_stack.top().first] = i - temper_stack.top().first;
            temper_stack.pop();
        }
        temper_stack.push({i,T[i]});
    }
    return results;
}

/*
给定状态下出现观测的概率
*/

std::vector<int> plusOne(std::vector<int>& digits) {
    std::vector<int> results;
    int over = 0;
    int num = 0;
    for(int i = (digits.size() - 1); i >=0 ; i--) {
        if(i == digits.size() - 1)
            num = digits.at(i) + 1 + over;
        else
        {
            num = digits.at(i) + over;
        }
        over = num / 10;
        results.push_back(num % 10);
    }
    if(over != 0)
        results.push_back(over);
    std::reverse(results.begin(),results.end());
    return results;    
}

int validIndex(const int i, const int j, const std::vector<std::vector<int>>& M) {
    //索引有效，返回像素值，否则为-1
    if(i >= 0 && i < M.size() && j >= 0 && j < M.front().size())
         return M[i][j];
    return -1;     
}

std::vector<std::vector<int>> imageSmoother(std::vector<std::vector<int>>& M) {
    int rows = M.size();
    int cols = M.front().size();
    std::vector<std::vector<int>> smooth_image(rows,std::vector<int>(cols,0));
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            int valid_count = 0;
            int sum = 0;
            for(int k = -1; k <=1; k++) {
                for(int l = -1; l <= 1; l++) {
                    if(validIndex(i+k,j+l,M) != -1) {
                        sum += validIndex(i+k,j+l,M);
                        valid_count++;
                    }
                }
            }
            std::cout<< sum <<std::endl;
            smooth_image[i][j] = std::floor(sum / valid_count);
        }
    }
    return smooth_image;
}

int majorityElement(std::vector<int>& nums) {
    int maj_num = 0;
    int maj_count = 0;
    std::unordered_map<int,int> has;
    for(auto& num : nums) {
        has[num] ++;
        if(has[num] > maj_count) {
            maj_count = has[num];
            maj_num = num;
        }
    }
    return maj_num;
}

struct ListNode {
    int val;
    ListNode* next;
    ListNode(int v) : val(v),next(nullptr) {}
};

// 分割链表，使得小于x的节点在大于等于x节点的左边
ListNode* partition(ListNode* head, int x) {
    ListNode* p = head;
    ListNode* dummysmall = nullptr;
    ListNode* dummylarge = nullptr;
    ListNode* dummyp = nullptr;
    ListNode* dummyq = nullptr;
    while(p != nullptr) {
        if(p->val < x) {
            if(dummysmall == nullptr) {
                dummysmall = p;
                dummyp = dummysmall;
            } else {
                dummyp->next = p;
                dummyp = p;
            }
        } else {
            if(dummylarge == nullptr) {
                dummylarge = p;
                dummyq = dummylarge;
            } else {
                dummyq->next = p;
                dummyq = p;
            }
        }
        //传递好之后，断开之前的节点
        p = p->next;
        if(dummyp)
            dummyp->next = nullptr;
        if(dummyq)
            dummyq->next = nullptr;    
    }
    if(dummysmall == nullptr)
        return dummylarge;
    dummyp->next = dummylarge;
    return dummysmall;    
}

//3000ms内的请求次数统计
class RecentCounter {
public:
    RecentCounter() {
        
    }
    
    int ping(int t) {
        int count = 1;
        for(int i = int(timestamps.size()-1); i >= 0; i--) {
            std::cout<<"i: "<<i<<std::endl;
            if(timestamps[i] >= (t - 3000))
                count++;
            else
                break;    
        }
        timestamps.push_back(t);
        return count;
    }

    std::vector<int> timestamps;
};

//基于线性表的循环队列实现
class MyCircularQueue {
public:
    /** Initialize your data structure here. Set the size of the queue to be k. */
    MyCircularQueue(int k) {
        data_nums = 0;
        size = k;
        front_index = 0;
        end_index = 0;
        queue_data = std::vector<int>(size,-1);
    }
    
    /** Insert an element into the circular queue. Return true if the operation is successful. */
    bool enQueue(int value) {
        if(isFull())
            return false;
        queue_data[end_index] = value;
        end_index = (end_index + 1) % size;
        data_nums++;
        return true;
    }
    
    /** Delete an element from the circular queue. Return true if the operation is successful. */
    bool deQueue() {
        if(isEmpty())
            return false;
        front_index = (front_index + 1) % size;
        data_nums--;
        return true;
    }
    
    /** Get the front item from the queue. */
    int Front() {
        if(isEmpty())
            return -1;
        return queue_data[front_index];
    }
    
    /** Get the last item from the queue. */
    int Rear() {
        if(isEmpty())
            return -1;
        int rear_index = end_index - 1;
        if(rear_index < 0)
            rear_index += size;
        return queue_data[rear_index];
    }
    
    /** Checks whether the circular queue is empty or not. */
    bool isEmpty() {
        return data_nums == 0;
    }
    
    /** Checks whether the circular queue is full or not. */
    bool isFull() {
        return data_nums == size;
    }

    std::vector<int> queue_data;
    int data_nums;
    int front_index;
    int end_index;
    int size;
};

// 循环双端队列
class MyCircularDeque {
public:
    /** Initialize your data structure here. Set the size of the deque to be k. */
    MyCircularDeque(int k) {
        queue_data = std::vector<int>(k,-1);
        data_nums = 0;
        front_index = -1;
        end_index = 0;
        size = k;
    }
    
    /** Adds an item at the front of Deque. Return true if the operation is successful. */
    bool insertFront(int value) {
        if(isFull())
            return false;
        if(front_index == -1) {
            if(isEmpty()) {
                queue_data[0] = value;
                front_index = 0; // front_index 指向当前的队首元素
                end_index = 1;
            } else {
                queue_data.back() = value;
                front_index = size - 1;
            }
        } else {
            front_index -= 1;
            if(front_index < 0)
                front_index += size;
            queue_data[front_index] = value;    
        }    
        data_nums++;
        return true;
    }
    
    /** Adds an item at the rear of Deque. Return true if the operation is successful. */
    bool insertLast(int value) {
        if(isFull())
            return false;
        queue_data[end_index] = value;
        end_index = (end_index + 1) % size;
        data_nums++;
        return true;
    }
    
    /** Deletes an item from the front of Deque. Return true if the operation is successful. */
    bool deleteFront() {
        if(isEmpty())
            return false;
        front_index = (front_index + 1) % size;
        data_nums--;
        return true; 
    }
    
    /** Deletes an item from the rear of Deque. Return true if the operation is successful. */
    bool deleteLast() {
        if(isEmpty())
            return false;
        end_index--;
        if(end_index < 0)
            end_index += size;
        data_nums--;
        return true;     
    }
    
    /** Get the front item from the deque. */
    int getFront() {
        if(isEmpty())
            return -1;
        return queue_data[front_index];    
    }
    
    /** Get the last item from the deque. */
    int getRear() {
        if(isEmpty())
            return -1;
        int index = end_index - 1;
        if(index < 0)
            index += size;
        return queue_data[index];        
    }
    
    /** Checks whether the circular deque is empty or not. */
    bool isEmpty() {
        return data_nums == 0;
    }
    
    /** Checks whether the circular deque is full or not. */
    bool isFull() {
        return data_nums == size;
    }

    std::vector<int> queue_data;
    int data_nums;
    int front_index;
    int end_index;
    int size;
};

// binary tree
// left == left_child right = right_child
std::vector<int> inorderTraversal(TreeNode* root) {
    if(!root)
        return {};
    std::vector<int> result = inorderTraversal(root->left_child);
    result.push_back(root->val);
    auto right_result = inorderTraversal(root->right_child);
    result.insert(result.end(),right_result.begin(),right_result.end());
    return result; 
}

bool isSameTree(TreeNode* p, TreeNode* q) {
    if(p == nullptr && q == nullptr)
        return true;
    if(p == nullptr & q != nullptr)
        return false;
    if(p != nullptr && q == nullptr)
        return false;
    if(p ->val != q->val)
        return false;
    return isSameTree(p->left_child,q->left_child) && 
           isSameTree(p->right_child,q->right_child);                        
}

//左树和右树是否形成对称结构
bool isSymmetric(TreeNode* left_tree, TreeNode* right_tree) {
    if(left_tree == nullptr && right_tree == nullptr)
        return true;
    if(left_tree == nullptr & right_tree != nullptr)
        return false;
    if(left_tree != nullptr && right_tree == nullptr)
        return false;
    if(left_tree->val != right_tree->val)
        return false;
    return isSymmetric(left_tree->left_child,right_tree->right_child) &&
           isSymmetric(left_tree->right_child,right_tree->left_child);        
}

bool isSymmetric(TreeNode* root) {
    if(root == nullptr)
        return true;
    return isSymmetric(root->left_child,root->right_child);    
}

//二叉树的最大深度
int maxDepth(TreeNode* root) {
    if(root == nullptr)
        return 0;
    return 1 + std::max(maxDepth(root->left_child), maxDepth(root->right_child));        
}

// 到叶子节点的最小深度
int minDepth(TreeNode* root) {
    if(root == nullptr)
        return 0;
    if(root->left_child == nullptr)
        return 1 + minDepth(root->right_child);
    if(root->right_child == nullptr)
        return 1 + minDepth(root->left_child);
    //左右都有才进行深度大小的比较            
return 1 + std::min(minDepth(root->left_child),minDepth(root->right_child));            
}

// 左右子树高度上是否平衡
// 先求每个节点上的深度
bool isBalanced(TreeNode* root) {
    if(root == nullptr)
        return true;
    int left_depth = maxDepth(root->left_child);
    int right_depth = maxDepth(root->right_child);
    if(std::abs(left_depth - right_depth) > 1)
        return false;
    return isBalanced(root->left_child) && isBalanced(root->right_child);    
}

// 最大平均数的连续子数组
double findMaxAverage(std::vector<int>& nums, int k) {
    // brute force sulation
    // int max_sum = std::numeric_limits<int>::min();
    // for(int i = 0; i <= nums.size() - k; i++) {
    //     int sum = 0;
    //     for(int j = i; j < i + k; j++) {
    //         sum += nums[j];
    //     }
    //     if(sum > max_sum) {
    //         max_sum = sum;
    //     }
    // }
    // return max_sum / double(k);

    // sliding window sulation
    int start_sum = 0;
    for(int i = 0; i < k; i++) {
        start_sum += nums[i];
    }
    int max_sum = start_sum;
    int start_window_index = 1;
    for(; start_window_index <= nums.size() - k; start_window_index++) {
        int tmp_sum = start_sum - nums[start_window_index - 1] + nums[start_window_index + k - 1];
        if(tmp_sum > max_sum)
            max_sum = tmp_sum;
        start_sum = tmp_sum; // 上一阶段的子数组之和    
    }
    return max_sum / double(k);
}

int trap(std::vector<int>& height) {
    std::vector<int> left_high(height.size(), 0);
    std::vector<int> right_high(height.size(),0);
    int left_max = 0;
    for(int i = 0; i < int(height.size()); i++) {
        left_high.at(i) = left_max;
        left_max = std::max(left_max,height[i]);
    }
    int right_max = 0;
    for(int  i = height.size() - 1; i >= 0; i--) {
        right_high[i] = right_max;
        right_max = std::max(right_max, height[i]);
    }

    int result = 0;
    for(int i = 1; i < height.size() - 1; i++) {
        int current_height = height[i];
        int left_max_height = left_high[i];
        int right_max_height = right_high[i];
        result += std::max(0, std::min(left_max_height,right_max_height) - current_height);
    }
    return result;
}

int findPeakElement(std::vector<int>& nums) {
    for(int i = 0; i < nums.size(); i++) {
        if(i == 0) {
            if(nums[i] > nums[i+1])
                return 0;
        }
        if(i == nums.size() - 1) {
            if(nums[i] > nums[i-1])
                return i-1;
        }
        if(nums[i] > nums[i-1] && nums[i] > nums[i+1])
            return i;
        continue;
    }
    return 0;
}

void zhushi() {
    #if 0
    std::cout<<"zhushi"<<std::endl;
    #endif

    #if 1
    std::cout<<"not zhushi"<<std::endl;
    #endif
}

// 联合体会覆盖内存空间值
union a_union {
    struct B {
        int x;
        int y;
    } b;

    int k;
};

enum weekday {
    sun, mon, the, wed, thu,may,sat
};

struct student
{
    int num;
    char name[20];
    char gender;
};

//使用类的前向声明定义时，只允许使用指针这样不完整的定义
class B;

class A {
public:

private:
B* b;
};

class Application
{ 
public:
    static void f(); 
    static void g();
private:
    static int global;
};

int Application::global=0;

void Application::f() {  global=5; }
void Application::g() {  std::cout<<global<<std::endl;}

// 常对象只能调用其对应的常成员函数
// void print() const; const A a;

// 根据实参调用来确定模版函数中参数的类型，然后编译器生成相应类型的模版函数
template<typename T>
void sort(T& a, int n){
    for(int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            if(a[j] > a[i]) {
                int tmp = a[i];
                a[i] = a[j];
                a[j] = tmp;
            }
        }
    }
}

template<typename T>
void display(T& a,int n){
    for(int i = 0; i < n;i++) {
        std::cout<<a[i]<<" ";
    }
    std::cout<<std::endl;
}

template<typename T, int MAXSIZE>
class Stack {
    private:
    int top = -1;
    T elem[MAXSIZE];

    public:
    bool empty() {
        if(top <= -1)
            return 1;
        return 0;    
    }

    bool full() {
        if(top >= (MAXSIZE-1))
            return 1;
        return 0;    
    }

    void push(T e);
    T pop();
};

// 模版类定义
template<typename T, int MAXSIZE>
void Stack<T, MAXSIZE>::push(T e) {
    if(full()) {
        std::cout<<"already full"<<std::endl;
        return;
    }
    elem[++top] = e;
}

template<typename T, int MAXSIZE>
T Stack<T, MAXSIZE>::pop() {
    if(empty()) {
        std::cout<<"already empty"<<std::endl;
        return;
    }
    return elem[top--];
}

template<typename T>
using Vec = std::vector<T, std::allocator<T>>;

//没有explicit声明，类对象会进行隐式转换
//加入explicit声明，类对象不会进行隐式转换 

//final 放在类后面表示该类不能被继承
//final 放在函数后面表示该函数不能被override

//virtual override会做重载虚函数检查

void findOdd(unsigned long long start, unsigned long long end) {
    unsigned long long odd_count = 0;
    for(unsigned long long i = start; i < end; i++) {
        if((i & 1) == 1)
            odd_count++;
    }
}

void findEven(unsigned long long start, unsigned long long end) {
    unsigned long long even_count = 0;
    for(unsigned long long i = start; i < end; i++) {
        if((i & 1) == 0)
            even_count++;
    }
}

void compare_time() {
    unsigned long long start = 0;
    unsigned long long end = 19000000;
    auto start_time = std::chrono::high_resolution_clock::now();
    findOdd(start, end);
    findEven(start, end);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = 
        std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout<<"un thread time: "<< duration.count() / 1000000.0 <<"s";

    start_time = std::chrono::high_resolution_clock::now();
    std::thread t1(findOdd, start, end);
    std::thread t2(findEven, start, end);
    end_time = std::chrono::high_resolution_clock::now();
    duration = 
        std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout<<"thread time: "<< duration.count() / 1000000.0 <<"s";
}

class Base {
    public:
    static void func(int n) {
        while(n--) {
            std::cout<< n << " ";
        }
    }
};

void run(int count) {
    while (count-- > 0) {
        std::cout << count << std::endl;
    }
    std::this_thread::sleep_for(std::chrono::seconds(3));
}

//数据竞争导致数据读写的错误
int sum = 0;
std::mutex m;
void countgold() {
    int i; //local to each thread
    for (i = 0; i < 10000000; i++) {
        sum += 1;
    }
}

int string2int(std::string str) {
    return std::stoi(str);
}

ListNode* reverseList(ListNode* head) {
    if(head == nullptr || head->next == nullptr)
        return head;
    ListNode* pre_node = head;
    ListNode* cur_node = pre_node->next;
    pre_node->next = nullptr;

    while(cur_node->next != nullptr) {
        ListNode* next_node = cur_node->next;
        cur_node->next = pre_node;
        pre_node = cur_node;
        cur_node = next_node;
    }
    cur_node->next = pre_node;
    return cur_node;    
}

ListNode* swapPairs(ListNode* head) {
    if(head == nullptr || head->next == nullptr)
        return head;
    ListNode* pre = head;
    ListNode* cur = pre->next;
    ListNode* result = head->next;
    ListNode* next = cur->next;

    while(pre != nullptr && cur != nullptr && next != nullptr) {
        cur->next = pre;
        if(next->next != nullptr)
            pre->next = next->next;
        else
        {
            pre->next = next;
        }
            

        pre = next;
        if(pre == nullptr)
            break;
        cur = pre->next;
        if(cur == nullptr)
            break;
        next = cur->next;
    }

    if(cur == nullptr)
        return result;

    if(next == nullptr) {
        cur->next = pre; 
        pre->next = nullptr;
    }
           

   return result;
        
}

bool hasCycle(ListNode *head) {
    if(head == nullptr || head->next == nullptr)
        return false;
    ListNode* Fast = head->next;
    ListNode* Slow = head;
    while(Fast != nullptr && Slow != nullptr) {
        if(Fast == Slow)
            return true;
        Slow = Slow->next;
        if(Fast->next == nullptr)
            return false;
        Fast = Fast->next->next;        
    }

    return false;
}

class KthLargest {
public:
    KthLargest(int k, std::vector<int>& nums) {
        k_ = k;
        for(int i : nums) {
            if(q.size() < k)
                q.push(i);
            else
            {
                if(i > q.top()) {
                    q.pop();
                    q.push(i);
                }
                    
            }
            
        }
    }
    
    int add(int val) {
        if(q.size() < k_) {
            q.push(val);
            return q.top();
        }
        if(val <= q.top())
            return q.top();
        q.pop();    
        q.push(val);
        return q.top();    
    }
    
    std::priority_queue<int, std::vector<int>, std::greater<int> > q;
    int k_;
};

/**
 * Your KthLargest object will be instantiated and called as such:
 * KthLargest* obj = new KthLargest(k, nums);
 * int param_1 = obj->add(val);
 */

//用双端队列deque求解sliding window问题
std::vector<int> maxSlidingWindow(std::vector<int>& nums, int k) {
    std::vector<int> results;
    return results;
}

bool isAnagram(std::string s, std::string t) {
    std::unordered_map<char, int> char_count;
    for(char c : s) {
        if(char_count.count(c) == 0)
            char_count[c] = 1;
        else
        {
            char_count[c] ++;
        }      
    }
    for(char c : t) {
        if(char_count[c] == 0)
            return false;
        else
        {
            char_count[c]--;
        }    
    }
    for(auto e : char_count) {
        if(e.second != 0)
            return false;
    }

    return true;
}

std::vector<std::vector<int>> threeSum_1(std::vector<int>& nums) {
    std::vector<std::vector<int>> results;

    std::map<int,int> num_count;
    for(int& num : nums) {
        if(num_count.count(num) == 0)
            num_count[num] = 1;
        else {
            num_count[num] += 1;
        }     
    }

    for(int i = 0; i < nums.size(); i++) {
        for(int j = i+1; j < nums.size(); j++) {
            int target_num = -(nums[i] + nums[j]);
            int count = 1;
            if(target_num == nums[i])
                count++;
            if(target_num == nums[j])
                count++;
            if(num_count[target_num] == count)
                results.push_back({nums[i], nums[j], target_num});        
        }
    }
    return results;
}

//中序遍历的二叉搜索树是一个升序的结构
void inorder_tranverce(TreeNode* root, std::vector<int>& results) {
    if(root == nullptr)
        return;
    inorder_tranverce(root->left_child, results);
    results.push_back(root->val);
    inorder_tranverce(root->right_child, results);            
}

bool isValidBST(TreeNode* root) {
    std::vector<int> inorder_results;
    inorder_tranverce(root, inorder_results);
    for(size_t i = 0; i < inorder_results.size() - 1; i++) {
        if(inorder_results[i+1] <= inorder_results[i])
            return false;
    }
    return true;
}

//二叉搜索树的最近公共祖先
TreeNode* lowestCommonAncestor_BST(TreeNode* root, TreeNode* p, TreeNode* q) {
    if(root == nullptr || p == nullptr || q == nullptr)
        return root;
    if(p->val < root->val && q->val > root->val)
        return root;
    if(p->val < root->val && q->val < root->val)
        return lowestCommonAncestor_BST(root->left_child, p, q);
    if(p->val > root->val && q->val > root->val)
        return lowestCommonAncestor_BST(root->right_child, p, q);
    return root;
        
}

//二叉树的最近公共祖先
TreeNode* lowestCommonAncestor(TreeNode* root, TreeNode* p, TreeNode* q) {
    if(root == nullptr || root == p || root == q)
        return root;
    TreeNode* left = lowestCommonAncestor(root->left_child, p, q);
    TreeNode* right = lowestCommonAncestor(root->right_child, p, q);
    if(left == nullptr)
        return right;
    if(right == nullptr)
        return left;
    return root;            
}

double myPow(double x, int n) {
    //return std::pow(x,n);
    if(n == 0 && x != 0)
        return 1;
    if(n < 0)
        return 1 / myPow(x, -n);
    double pow = 1;
    if(n % 2 == 1)
        return x * myPow(x * x, n / 2);
    return myPow(x * x, n / 2);      
}

//二叉树层次遍历输出
std::vector<std::vector<int>> levelOrder(TreeNode* root) {
    std::vector<std::vector<int>> results;
    if(root == nullptr)
        return results;
    std::queue<TreeNode*> q;
    q.push(root);
    while(!q.empty()) {
        std::vector<int> result;
        std::vector<TreeNode*> current_level_nodes;
        //记录当前层级的节点
        while(!q.empty()) {
            result.push_back(q.front()->val);
            current_level_nodes.push_back(q.front());
            q.pop();
        }
        results.push_back(result);
        for(auto& node : current_level_nodes) {
            if(node->left_child != nullptr)
                q.push(node->left_child);
            if(node->right_child != nullptr)
                q.push(node->right_child);    
        }
    }
    return results;    
}

//left right代表已经用了的括号的个数
void gen_helper(int n, int left, int right, std::string str ,std::vector<std::string>& results) {
    if(left == n && right == n) {
        results.push_back(str);
        return;
    }

    if(left < n) {
        gen_helper(n, left+1, right, str + "(", results);
    }
    if(right < n && right < left) {
        gen_helper(n, left, right+1,str + ")",results);
    }
}

//生成有效的括号集合
std::vector<std::string> generateParenthesis(int n) {
    std::vector<std::string> results;
    std::string str = "";
    gen_helper(n, 0, 0, str, results);
    return results;
}

int main()
{
    std::vector<std::string> results = generateParenthesis(3);
    for(auto& result : results) {
        std::cout<< result << std::endl;
    }
    return 0;
}
