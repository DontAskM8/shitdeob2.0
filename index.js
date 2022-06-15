const fs = require("fs")
const beautify = require("js-beautify")
const parser = require("@babel/parser")
const traverse = require("@babel/traverse").default
const generator = require("@babel/generator").default
const types = require("@babel/types")

var original = fs.readFileSync("./raw.js", "utf8")

var undefinedShit = /''\+.*?(?=\))/g
for(let i of original.match(undefinedShit)){
	original = original.replace(i, "undefined")
}


// var evaluableMatches = original.match(/[A-z]+=.*?(,|;)/g)
var evaluableMatches = /(([A-z]|[0-9])+)=(((((\[|\]|\+|!|-| )+)(?=(,|;))))|(((([A-z]|[0-9])+(\*|\+|-|\/)+)+([A-z]|[0-9])+)(?=(,|;))))/g

function blEval(x) {
	var arr = ["{}", "[]", "this", "global"]
	for(let i of arr){
		if((x == i) || (x == (i + ",")) || (x == (i + ";"))){
			return true
		}
	}
	return false
}


// Until line 65 shall be updated to babel but im lazy
var functionsToExecute = /(?<!\w)([A-z]|[0-9])+\(\)/g

for(let k of original.match(functionsToExecute)){
	var contentOfFunction = new RegExp(`(?<=function ${k}\\(\\){).*?(?=})`, "g")
	for(let j of original?.match(contentOfFunction) ?? ""){
		for(let i of j.match(evaluableMatches) ?? ""){
			try {
				var evaled = eval(i)
				if(!blEval(i) && (evaled.toString() != "[object Object]") && !Array.isArray(evaled)){
					eval("global." + i)
					original = original.replace(i, [i.split("=")[0], "=", typeof evaled == "string" ? `'${evaled}'` : evaled].join(""))
				}
			}catch(e){
			}
		}
		
	}
}


for(let i of original.match(evaluableMatches) ?? ""){
	try {
		var evaled = eval(i)
		if(!blEval(i) && (evaled.toString() != "[object Object]") && !Array.isArray(evaled)){
			eval("global." + i)
			original = original.replace(i, [i.split("=")[0], "=", typeof evaled == "string" ? `'${evaled}'` : evaled].join(""))
		}
	}catch(e){
	}
}

var globalVars = Object.keys(global).filter(x => !['global', 'clearInterval', 'clearTimeout', 'setInterval', 'setTimeout', 'queueMicrotask', 'performance', 'clearImmediate', 'setImmediate'].includes(x))
function isGlobalVars(va){
	return globalVars.includes(va) && typeof global[va] != "function"
}

//Main function that triggers the entire universe
var executingFunction = /(?<=(return ))([A-z]|[0-9])+\.call\(this,([A-z]|[0-9])+\)/g
var [callingFunc, callingVal] = original.match(executingFunction)[0].replace(/call\(this,|\)/g, "").split(".")

console.log(callingFunc, callingVal)


var _testCode = parser.parse(original)
traverse(_testCode, {
	// Replace switch cases var to value
	SwitchStatement: function(path){
		var { node } = path
		node.cases = node.cases.map(x => {
			if(!isGlobalVars(x.test?.name)) return x;
			else {
				x.test = types.numericLiteral(global[x.test.name])
				return x
			}
		})
	},
	//Replace var refer var to value
	AssignmentExpression: function(path){
		var { left, right } = path.node
		if(left.type != "Identifier" || right.type != "Identifier" ) return;
		if(!globalVars.includes(right.name)) return;
		try {
			Object.assign(right, types.numericLiteral(global[right.name]))
		}
		catch{}
	},
	//Eval all functions make it global
	"FunctionDeclaration|VariableDeclaration": function(path){
		try {
			global[path.node.id.name] = eval(original.slice(path.node.start, path.node.end))
		}catch{
		}
		
		if(types.isVariableDeclaration(path.node) && path.node.declarations.length == 1){
			let { id, init } = path.node.declarations[0]
			if(id?.name != callingFunc) return;
			
			let callingParamFunc = init.id.name
			let callingParams = init.params.map(x => x.name)
			
			let whileLoopInMain = init.body.body.find(x => x.type == "WhileStatement")
			let whileCond = original.slice(whileLoopInMain.test.loc.start.index, whileLoopInMain.test.loc.end.index)
			
			
			
			function loopFunction(startVal){
				eval(`var ${callingParams[0]} = ${startVal}`)
				while(eval(whileCond)){
					// console.log(whileLoopInMain.body.body[0])
					let cases = whileLoopInMain.body.body[0].cases
					// console.log(cases[0])
					cases = cases.find(x => {
						if(types.isIdentifier(x.test)){
							return eval(x.test?.name) == eval(callingParams[0])
						}else{
							console.log(x.test)
							return x.test?.value == eval(callingParams[0])
						}
					})
					
					let { consequent } = cases
					consequent = consequent.find(x => types.isBlockStatement(x))
					
					let { body } = consequent
					body = body.filter(x => types.isExpressionStatement(x))
					for(let content of body){
						let { expression } = content	
						//Execute the same loop
						if(types.isCallExpression(expression) && expression.callee.name == callingParamFunc){
							loopFunction(eval(expression.arguments.find(x => types.isIdentifier(x)).name))
							
						}
						else if(types.isAssignmentExpression(expression)){
							if(expression.left.name == callingParams[0]){
								eval(original.slice(expression.start, expression.end))
								console.log(`Loop ${eval(startVal)}: ${eval(callingParams[0])}`)
							}else if(expression.operator == "="){
								try {
									eval("global." + original.slice(expression.start, expression.end))
									if(typeof global[expression.left.name] == "number"){
										expression.right = types.NumericLiteral(global[expression.left.name])
									}
								}
								catch(e){
								}
							}else{
								try {
									eval(original.slice(expression.start, expression.end))
								}
								catch{}
							}
							
						}
					}
					if(!eval(whileCond)){
						console.log("Exiting loop", eval(startVal))
					}
				}
			}
			loopFunction(callingVal)
		}
	},
	SwitchCase: function(path){
		var { container, parent } = path
		var subject = parent.discriminant.name
		for(let cases of container){
			var testVal = cases.test?.value
			// console.log(testVal ?? original.slice(cases.loc.start.index, cases.loc.end.index))
			let block = cases.consequent.find(x => x.type == "BlockStatement")
			
		}
	},
})

original = generator(_testCode, {
	minified: true
}, original).code


var beautifiedCode = beautify(original, {
	indent_size: "\t",
	space_in_empty_paren: true,
	unescape_strings: true
})

fs.writeFileSync("./deobed.js", original)
fs.writeFileSync("./beautified.js", beautifiedCode)